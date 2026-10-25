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

#include "cov/topology_graph.hpp"
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>

namespace cov {
namespace {
using Adjacency=std::vector<std::vector<std::uint32_t>>;
Adjacency graph(std::uint32_t count,const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges) {
    Adjacency result(count);
    for(const auto& [a,b]:edges) {
        if(a>=count || b>=count || a==b)throw std::invalid_argument("Invalid simple graph edge");
        result[a].push_back(b);result[b].push_back(a);
    }
    for(auto& row:result) {
        std::sort(row.begin(),row.end());row.erase(std::unique(row.begin(),row.end()),row.end());
    }
    return result;
}
struct Edge { std::size_t to,reverse; int capacity; bool original; };
using Network=std::vector<std::vector<Edge>>;
void add(Network& network,std::size_t a,std::size_t b) {
    const auto back_a=network[a].size(),back_b=network[b].size();
    network[a].push_back({b,back_b,1,true});network[b].push_back({a,back_a,0,false});
}
std::optional<std::array<std::vector<std::uint32_t>,2>> paths(
    const Adjacency& adjacency,std::uint32_t hub,
    const std::array<std::uint32_t,2>& sources,const std::array<std::uint32_t,2>& sinks) {
    const auto n=adjacency.size(),source=2*n,sink=source+1;
    Network network(2*n+2);
    for(std::size_t v=0;v<n;++v) {
        if(v==hub)continue;
        add(network,2*v,2*v+1);
        for(const auto neighbour:adjacency[v])if(neighbour!=hub)add(network,2*v+1,2*neighbour);
    }
    for(const auto v:sources)add(network,source,2*v);
    for(const auto v:sinks)add(network,2*v+1,sink);
    const auto absent=std::numeric_limits<std::size_t>::max();
    auto route=[&](bool positive_flow) {
        std::vector<std::pair<std::size_t,std::size_t>> parent(network.size(),{absent,absent});
        std::queue<std::size_t> queue;queue.push(source);parent[source]={source,0};
        while(!queue.empty() && parent[sink].first==absent) {
            const auto a=queue.front();queue.pop();
            for(std::size_t i=0;i<network[a].size();++i) {
                const auto& edge=network[a][i];
                const bool usable=positive_flow?(edge.original && edge.capacity==0):edge.capacity>0;
                if(usable && parent[edge.to].first==absent) {parent[edge.to]={a,i};queue.push(edge.to);}
            }
        }
        return parent;
    };
    // Unit vertex capacities require only two integral augmentations.
    for(int unit=0;unit<2;++unit) {
        const auto parent=route(false);if(parent[sink].first==absent)return std::nullopt;
        for(auto b=sink;b!=source;) {
            const auto [a,i]=parent[b];auto& edge=network[a][i];
            --edge.capacity;++network[b][edge.reverse].capacity;b=a;
        }
    }
    std::array<std::vector<std::uint32_t>,2> result;
    for(auto& path:result) {
        const auto parent=route(true);
        if(parent[sink].first==absent)throw std::logic_error("Integral flow decomposition failed");
        for(auto b=sink;b!=source;) {
            if(b<2*n && b%2==0)path.push_back(static_cast<std::uint32_t>(b/2));
            const auto [a,i]=parent[b];++network[a][i].capacity;b=a;
        }
        std::reverse(path.begin(),path.end());
    }
    return result;
}
}

bool verify_cycle_union(std::uint32_t count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,const CycleUnionWitness& witness) {
    const auto adjacency=graph(count,edges);if(witness.hub>=count)return false;
    std::array<std::set<std::uint32_t>,2> vertices;
    for(std::size_t j=0;j<2;++j) {
        const auto& cycle=witness.cycles[j];
        if(cycle.size()<4 || cycle.front()!=witness.hub || cycle.back()!=witness.hub)return false;
        for(std::size_t i=0;i+1<cycle.size();++i) {
            const auto a=cycle[i],b=cycle[i+1];
            if(a>=count || b>=count || !vertices[j].insert(a).second ||
                !std::binary_search(adjacency[a].begin(),adjacency[a].end(),b))return false;
        }
    }
    for(const auto v:vertices[0])if(v!=witness.hub && vertices[1].count(v))return false;
    return true;
}

std::optional<CycleUnionWitness> find_cycle_union_at_hub(std::uint32_t count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,std::uint32_t hub) {
    if(hub>=count)throw std::invalid_argument("Invalid cycle-union hub");
    const auto adjacency=graph(count,edges);const auto& neighbours=adjacency[hub];
    for(std::size_t a=0;a<neighbours.size();++a)
    for(std::size_t b=a+1;b<neighbours.size();++b)
    for(std::size_t c=b+1;c<neighbours.size();++c)
    for(std::size_t d=c+1;d<neighbours.size();++d) {
        const std::array<std::uint32_t,4> ends{{neighbours[a],neighbours[b],neighbours[c],neighbours[d]}};
        for(std::size_t other=1;other<4;++other) {
            std::array<std::uint32_t,2> sources{{ends[0],ends[other]}},sinks{};std::size_t index=0;
            for(std::size_t i=1;i<4;++i)if(i!=other)sinks[index++]=ends[i];
            const auto found=paths(adjacency,hub,sources,sinks);if(!found)continue;
            CycleUnionWitness witness;witness.hub=hub;
            for(std::size_t i=0;i<2;++i) {
                auto& cycle=witness.cycles[i];cycle.push_back(hub);
                cycle.insert(cycle.end(),(*found)[i].begin(),(*found)[i].end());cycle.push_back(hub);
            }
            if(!verify_cycle_union(count,edges,witness))throw std::logic_error("Constructed cycles are not a valid union");
            std::vector<bool> reached(count,false);std::queue<std::uint32_t> queue;
            for(const auto v:(*found)[0]) {reached[v]=true;queue.push(v);}
            while(!queue.empty()) {
                const auto v=queue.front();queue.pop();
                for(const auto next:adjacency[v])if(next!=hub && !reached[next]) {reached[next]=true;queue.push(next);}
            }
            for(const auto v:(*found)[1])witness.additional_connection_without_hub|=reached[v];
            return witness;
        }
    }
    return std::nullopt;
}

ChannelCycleSearchResult find_channel_cycle_union_at_hub(std::uint32_t count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,
    std::uint32_t hub,const std::array<std::vector<std::uint32_t>,2>& channels) {
    if(hub>=count)throw std::invalid_argument("Invalid channel-cycle hub");
    const auto adjacency=graph(count,edges);
    std::array<std::vector<bool>,2> member{{std::vector<bool>(count),std::vector<bool>(count)}};
    for(std::size_t j=0;j<2;++j)for(const auto atom:channels[j]) {
        if(atom>=count)throw std::invalid_argument("Invalid channel atom");
        member[j][atom]=true;
    }
    ChannelCycleSearchResult result;
    for(std::size_t j=0;j<2;++j)for(std::uint32_t atom=0;atom<count;++atom)
        if(atom!=hub && member[j][atom] && !member[1-j][atom])++result.exclusive_channel_atoms[j];
    const auto unscoped=find_cycle_union_at_hub(count,edges,hub);
    result.unscoped_cycle_union_available=unscoped.has_value();
    if(!unscoped || result.exclusive_channel_atoms[0]<2 || result.exclusive_channel_atoms[1]<2)return result;

    // Enumerate only cycles consistent with the declared channel boundary.
    // The preceding polynomial graph test rejects impossible unions first.
    // There is no finite cycle-count cutoff that could silently mean "absent".
    using Cycle=std::vector<std::uint32_t>;
    auto find_cycle=[&](std::size_t channel,const std::vector<bool>& forbidden,
                        const std::function<bool(const Cycle&)>& accept) {
        std::vector<bool> used=forbidden;
        Cycle path{hub};used[hub]=true;
        std::function<bool(std::uint32_t,std::size_t)> visit;
        visit=[&](std::uint32_t current,std::size_t matched) {
            for(const auto next:adjacency[current]) {
                if(next==hub) {
                    if(path.size()>=3 && path[1]<path.back() && matched>=2) {
                        Cycle closed=path;closed.push_back(hub);
                        if(accept(closed))return true;
                    }
                } else if(!used[next] && !member[1-channel][next]) {
                    used[next]=true;path.push_back(next);
                    const bool found=visit(next,matched+(member[channel][next]?1u:0u));
                    path.pop_back();used[next]=false;
                    if(found)return true;
                }
            }
            return false;
        };
        return visit(hub,0);
    };
    find_cycle(0,std::vector<bool>(count),[&](const Cycle& first) {
        ++result.qualifying_first_cycles_examined;
        std::vector<bool> forbidden(count);
        for(const auto atom:first)if(atom!=hub)forbidden[atom]=true;
        return find_cycle(1,forbidden,[&](const Cycle& second) {
            CycleUnionWitness witness;witness.hub=hub;witness.cycles={{first,second}};
            if(!verify_cycle_union(count,edges,witness))throw std::logic_error("Invalid channel-scoped cycle union");
            std::vector<bool> reached(count);std::queue<std::uint32_t> queue;
            for(const auto atom:first)if(atom!=hub) {reached[atom]=true;queue.push(atom);}
            while(!queue.empty()) {
                const auto atom=queue.front();queue.pop();
                for(const auto next:adjacency[atom])if(next!=hub && !reached[next]) {
                    reached[next]=true;queue.push(next);
                }
            }
            for(const auto atom:second)if(atom!=hub)witness.additional_connection_without_hub|=reached[atom];
            result.witness=std::move(witness);return true;
        });
    });
    return result;
}
}
