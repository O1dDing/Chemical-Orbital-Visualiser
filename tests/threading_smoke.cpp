/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#include "cov/threading.hpp"

#include <atomic>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace {

void require(bool condition, const char* name) {
    if (!condition) throw std::runtime_error(name);
}

void stop_source_round_trip() {
    cov::stop_source source;
    const auto token=source.get_token();
    require(!token.stop_requested(),"token starts running");
    require(source.request_stop(),"first stop request changes state");
    require(token.stop_requested(),"token observes stop");
    require(!source.request_stop(),"second stop request is unchanged");
}

void destructor_requests_stop_and_joins() {
    std::atomic<bool> finished{false};
    {
        cov::jthread worker([&](cov::stop_token token) {
            while (!token.stop_requested()) std::this_thread::yield();
            finished.store(true,std::memory_order_release);
        });
        require(worker.joinable(),"worker starts joinable");
    }
    require(finished.load(std::memory_order_acquire),"destructor joined stopped worker");
}

void move_assignment_stops_previous_worker() {
    std::atomic<bool> old_finished{false},new_finished{false};
    cov::jthread old([&](cov::stop_token token) {
        while (!token.stop_requested()) std::this_thread::yield();
        old_finished.store(true,std::memory_order_release);
    });
    cov::jthread next([&](cov::stop_token token) {
        while (!token.stop_requested()) std::this_thread::yield();
        new_finished.store(true,std::memory_order_release);
    });
    old=std::move(next);
    require(old_finished.load(std::memory_order_acquire),
            "move assignment stopped and joined old worker");
    require(!next.joinable() && old.joinable(),"move assignment transferred worker");
    require(old.request_stop(),"moved worker accepts stop");
    old.join();
    require(new_finished.load(std::memory_order_acquire),"moved worker joined");
}

void ordinary_callable_and_move_constructor() {
    std::atomic<int> result{0};
    cov::jthread first([&](int value) { result.store(value,std::memory_order_release); },42);
    cov::jthread second(std::move(first));
    require(!first.joinable() && second.joinable(),"move constructor transferred worker");
    second.join();
    require(result.load(std::memory_order_acquire)==42,"callable without stop token");
}

} // namespace

int main() {
    try {
        stop_source_round_trip();
        destructor_requests_stop_and_joins();
        move_assignment_stops_previous_worker();
        ordinary_callable_and_move_constructor();
        std::cout << "threading_smoke: stop, join and move verified\n";
    } catch (const std::exception& error) {
        std::cerr << "threading_smoke: " << error.what() << '\n';
        return 1;
    }
}
