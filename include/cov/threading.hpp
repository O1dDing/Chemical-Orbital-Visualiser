#pragma once

#include <thread>

#if defined(__cpp_lib_jthread) && __cpp_lib_jthread >= 201911L && \
    !defined(COV_FORCE_PORTABLE_THREADS)
#include <stop_token>

namespace cov {
using stop_token = std::stop_token;
using stop_source = std::stop_source;
using jthread = std::jthread;
}

#else

#include <atomic>
#include <functional>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

namespace cov {

class stop_token {
public:
    stop_token() noexcept = default;
    [[nodiscard]] bool stop_requested() const noexcept {
        return flag_ && flag_->load(std::memory_order_acquire);
    }
    [[nodiscard]] bool stop_possible() const noexcept { return static_cast<bool>(flag_); }

private:
    explicit stop_token(std::shared_ptr<std::atomic<bool>> flag) noexcept
        : flag_(std::move(flag)) {}
    std::shared_ptr<std::atomic<bool>> flag_;
    friend class stop_source;
};

class stop_source {
public:
    stop_source() : flag_(std::make_shared<std::atomic<bool>>(false)) {}
    [[nodiscard]] stop_token get_token() const noexcept { return stop_token(flag_); }
    bool request_stop() noexcept {
        return flag_ && !flag_->exchange(true,std::memory_order_acq_rel);
    }
    [[nodiscard]] bool stop_requested() const noexcept {
        return flag_ && flag_->load(std::memory_order_acquire);
    }

private:
    std::shared_ptr<std::atomic<bool>> flag_;
};

class jthread {
public:
    jthread() = default;

    template <class F, class... Args,
              std::enable_if_t<!std::is_same_v<std::decay_t<F>,jthread>,int> = 0>
    explicit jthread(F&& function, Args&&... arguments)
        : source_(), worker_([token=source_.get_token(),
                              callable=std::decay_t<F>(std::forward<F>(function)),
                              values=std::make_tuple(std::forward<Args>(arguments)...)]() mutable {
            std::apply([&](auto&&... unpacked) {
                if constexpr (std::is_invocable_v<decltype(callable),stop_token,
                                                  decltype(unpacked)...>) {
                    std::invoke(std::move(callable),token,
                                std::forward<decltype(unpacked)>(unpacked)...);
                } else {
                    std::invoke(std::move(callable),
                                std::forward<decltype(unpacked)>(unpacked)...);
                }
            },std::move(values));
        }) {}

    ~jthread() { stop_and_join(); }
    jthread(const jthread&) = delete;
    jthread& operator=(const jthread&) = delete;
    jthread(jthread&& other) noexcept
        : source_(std::move(other.source_)), worker_(std::move(other.worker_)) {}
    jthread& operator=(jthread&& other) noexcept {
        if (this != &other) {
            stop_and_join();
            source_ = std::move(other.source_);
            worker_ = std::move(other.worker_);
        }
        return *this;
    }

    [[nodiscard]] bool joinable() const noexcept { return worker_.joinable(); }
    void join() { worker_.join(); }
    bool request_stop() noexcept { return source_.request_stop(); }
    [[nodiscard]] stop_token get_stop_token() const noexcept { return source_.get_token(); }

private:
    void stop_and_join() noexcept {
        if (worker_.joinable()) {
            source_.request_stop();
            worker_.join();
        }
    }
    stop_source source_;
    std::thread worker_;
};

} // namespace cov

#endif
