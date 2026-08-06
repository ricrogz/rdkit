/*
 * This is an extensive modification by Greg Landrum of
 * pieces from several files in the vflib-2.0 distribution
 *
 * The initial version of the modifications was completed
 *   in April 2009.
 *
 * the original author of the vflib files is:
 *    Author: P. Foggia
 *  http://amalfi.dis.unina.it/graph/db/vflib-2.0/doc/vflib.html
 *
 */
#include <RDGeneral/ControlCHandler.h>

#include <boost/graph/adjacency_list.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <queue>
#include <vector>

#ifndef __BGL_VF2_SUB_STATE_H__
#define __BGL_VF2_SUB_STATE_H__
#define RDK_ADJ_ITER typename Graph::adjacency_iterator

namespace boost {
namespace detail {
typedef std::uint32_t node_id;
const node_id NULL_NODE = 0xFFFFFFFF;
struct NodeInfo {
  node_id id;
  node_id in;
  node_id out;
  node_id candidates;
};

template <class Graph>
struct Pair {
  node_id n1, n2;
  bool hasiter{false};
  RDK_ADJ_ITER nbrbeg, nbrend;

  Pair() : n1(NULL_NODE), n2(NULL_NODE) {}
};

/**
 * The ordering by in/out degree
 */
static bool nodeInfoComp1(const NodeInfo &a, const NodeInfo &b) {
  if (a.out < b.out) {
    return true;
  }
  if (a.out > b.out) {
    return false;
  }
  if (a.in < b.in) {
    return true;
  }
  if (a.in > b.in) {
    return false;
  }
  return a.id < b.id;
}

/**
 * The ordering by frequency/valence.
 * The frequency is in the out field, the valence in `in'.
 */
static bool nodeInfoComp2(const NodeInfo &a, const NodeInfo &b) {
  if (!a.in && b.in) {
    return false;
  }
  if (a.in && !b.in) {
    return true;
  }
  if (a.candidates < b.candidates) {
    return true;
  }
  if (a.candidates > b.candidates) {
    return false;
  }
  if (a.out < b.out) {
    return true;
  }
  if (a.out > b.out) {
    return false;
  }
  if (a.in < b.in) {
    return true;
  }
  if (a.in > b.in) {
    return false;
  }
  return a.id < b.id;
}

template <class Graph, class VertexDescr, class EdgeDescr>
VertexDescr getOtherIdx(const Graph &g, const EdgeDescr &edge,
                        const VertexDescr &vertex) {
  VertexDescr tmp = boost::source(edge, g);
  if (tmp == vertex) {
    tmp = boost::target(edge, g);
  }
  return tmp;
}

struct CandidateMatrix {
  std::size_t targetSize{0};
  std::vector<unsigned char> values;
  std::vector<unsigned int> rowCounts;

  explicit operator bool() const { return !values.empty(); }

  bool get(node_id query, node_id target) const {
    return values[static_cast<std::size_t>(query) * targetSize + target] != 0;
  }

  void remove(node_id query, node_id target) {
    const auto idx = static_cast<std::size_t>(query) * targetSize + target;
    if (values[idx]) {
      values[idx] = 0;
      --rowCounts[query];
    }
  }

  bool hasEmptyRow() const {
    return std::find(rowCounts.begin(), rowCounts.end(), 0) != rowCounts.end();
  }
};

enum class NeighborhoodMatchResult {
  Match,
  NoMatch,
  Incomplete
};

template <class Edge>
struct IncidentEdge {
  node_id neighbor;
  Edge edge;
};

struct NeighborhoodWorkspace {
  std::vector<unsigned char> compatible;
  std::vector<unsigned int> queryOrder;
  std::vector<unsigned int> compatibleCounts;
  std::vector<int> targetMatches;
  std::vector<unsigned char> seen;
};

template <class Graph>
std::vector<std::vector<IncidentEdge<typename Graph::edge_descriptor>>>
BuildIncidentEdgeCache(const Graph &graph) {
  using Edge = typename Graph::edge_descriptor;
  const auto graphSize = boost::num_vertices(graph);
  std::vector<std::vector<IncidentEdge<Edge>>> result(graphSize);
  for (node_id vertex = 0; vertex < graphSize; ++vertex) {
    auto &incidentEdges = result[vertex];
    incidentEdges.reserve(boost::out_degree(vertex, graph));
    typename Graph::out_edge_iterator edgeBeg, edgeEnd;
    boost::tie(edgeBeg, edgeEnd) = boost::out_edges(vertex, graph);
    while (edgeBeg != edgeEnd) {
      incidentEdges.push_back(
          {static_cast<node_id>(getOtherIdx(graph, *edgeBeg, vertex)),
           *edgeBeg});
      ++edgeBeg;
    }
  }
  return result;
}

template <class Edge, class EdgeCompatible>
NeighborhoodMatchResult HasCompatibleNeighborhood(
    const std::vector<IncidentEdge<Edge>> &queryEdges,
    const std::vector<IncidentEdge<Edge>> &targetEdges,
    const CandidateMatrix &candidateMatrix, EdgeCompatible &edgeCompatible,
    std::size_t &work, std::size_t maxWork, NeighborhoodWorkspace &workspace) {
  if (queryEdges.size() > targetEdges.size()) {
    return NeighborhoodMatchResult::NoMatch;
  }

  const auto numTargetEdges = targetEdges.size();
  // Bound both the relation work and its temporary allocation. This
  // optimization is intended for the low-degree, simple graphs used by
  // molecular matching; high-degree generic graphs fall back to normal VF2.
  if (numTargetEdges && queryEdges.size() > (maxWork - work) / numTargetEdges) {
    return NeighborhoodMatchResult::Incomplete;
  }
  const auto relationSize = queryEdges.size() * numTargetEdges;
  work += relationSize;

  workspace.compatible.assign(relationSize, 0);
  workspace.queryOrder.resize(queryEdges.size());
  workspace.compatibleCounts.assign(queryEdges.size(), 0);
  auto &compatible = workspace.compatible;
  auto &queryOrder = workspace.queryOrder;
  auto &compatibleCounts = workspace.compatibleCounts;
  for (unsigned int queryEdgeIdx = 0; queryEdgeIdx < queryEdges.size();
       ++queryEdgeIdx) {
    queryOrder[queryEdgeIdx] = queryEdgeIdx;
    const auto queryNbr = queryEdges[queryEdgeIdx].neighbor;
    for (unsigned int targetEdgeIdx = 0; targetEdgeIdx < targetEdges.size();
         ++targetEdgeIdx) {
      const auto targetNbr = targetEdges[targetEdgeIdx].neighbor;
      if (candidateMatrix.get(queryNbr, targetNbr) &&
          edgeCompatible(queryEdges[queryEdgeIdx].edge,
                         targetEdges[targetEdgeIdx].edge)) {
        compatible[static_cast<std::size_t>(queryEdgeIdx) * numTargetEdges +
                   targetEdgeIdx] = 1;
        ++compatibleCounts[queryEdgeIdx];
      }
    }
    if (!compatibleCounts[queryEdgeIdx]) {
      return NeighborhoodMatchResult::NoMatch;
    }
  }

  std::sort(queryOrder.begin(), queryOrder.end(),
            [&compatibleCounts](unsigned int lhs, unsigned int rhs) {
              return compatibleCounts[lhs] < compatibleCounts[rhs];
            });

  workspace.targetMatches.assign(numTargetEdges, -1);
  auto &targetMatches = workspace.targetMatches;
  auto augment = [&](auto &self, unsigned int queryEdgeIdx,
                     std::vector<unsigned char> &seen) -> bool {
    for (unsigned int targetEdgeIdx = 0; targetEdgeIdx < numTargetEdges;
         ++targetEdgeIdx) {
      const auto compatIdx =
          static_cast<std::size_t>(queryEdgeIdx) * numTargetEdges +
          targetEdgeIdx;
      if (!compatible[compatIdx] || seen[targetEdgeIdx]) {
        continue;
      }
      seen[targetEdgeIdx] = 1;
      if (targetMatches[targetEdgeIdx] < 0 ||
          self(self, static_cast<unsigned int>(targetMatches[targetEdgeIdx]),
               seen)) {
        targetMatches[targetEdgeIdx] = static_cast<int>(queryEdgeIdx);
        return true;
      }
    }
    return false;
  };

  workspace.seen.resize(numTargetEdges);
  auto &seen = workspace.seen;
  for (const auto queryEdgeIdx : queryOrder) {
    std::fill(seen.begin(), seen.end(), 0);
    if (!augment(augment, queryEdgeIdx, seen)) {
      return NeighborhoodMatchResult::NoMatch;
    }
  }
  return NeighborhoodMatchResult::Match;
}

template <class Graph, class VertexCompatible, class EdgeCompatible>
CandidateMatrix BuildCandidateMatrix(const Graph &g1, const Graph &g2,
                                     VertexCompatible &vertexCompatible,
                                     EdgeCompatible &edgeCompatible,
                                     bool exactMatch) {
  constexpr std::size_t minQuerySize = 32;
  constexpr std::size_t maxCandidateCells = 1U << 20;
  constexpr std::size_t absoluteMaxRefinementWork = 8U << 20;
  constexpr std::size_t maxCachedEdges = 1U << 18;
  constexpr unsigned int maxRefinementPasses = 32;
  constexpr unsigned int maxUnrefinedRootCandidates = 8;

  CandidateMatrix result;
  const auto querySize = boost::num_vertices(g1);
  const auto targetSize = boost::num_vertices(g2);
  // The eager relation is reserved for large graph-isomorphism searches. It
  // solves the highly symmetric equal-size case without adding quadratic
  // setup to ordinary substructure searches (including repeated FMCS calls).
  if (!exactMatch || querySize < minQuerySize || !targetSize ||
      querySize > maxCandidateCells / targetSize) {
    return result;
  }

  result.targetSize = targetSize;
  result.values.resize(querySize * targetSize, 0);
  result.rowCounts.resize(querySize, 0);
  // Compatibility functors are BinaryPredicates: their answers must remain
  // stable while the relation is built and reused by both VF2 traversals.
  using Degree = typename boost::graph_traits<Graph>::degree_size_type;
  std::vector<Degree> targetDegrees(targetSize);
  for (node_id target = 0; target < targetSize; ++target) {
    targetDegrees[target] = boost::out_degree(target, g2);
  }

  std::size_t compatibilityChecks = 0;
  for (node_id query = 0; query < querySize; ++query) {
    const auto queryDegree = boost::out_degree(query, g1);
    for (node_id target = 0; target < targetSize; ++target) {
      if (!(compatibilityChecks++ & 0x3ffU) &&
          RDKit::ControlCHandler::getGotSignal()) {
        return result;
      }
      if (queryDegree == targetDegrees[target] &&
          vertexCompatible(query, target)) {
        result.values[static_cast<std::size_t>(query) * targetSize + target] =
            1;
        ++result.rowCounts[query];
      }
    }
    if (!result.rowCounts[query]) {
      return result;
    }
    // The legacy traversal starts at query vertex zero. If that root already
    // has few choices, normal VF2 avoids the pathological branching and is
    // cheaper than completing a quadratic compatibility relation.
    if (!query && result.rowCounts[query] <= maxUnrefinedRootCandidates) {
      return CandidateMatrix{};
    }
  }

  const auto maxRefinementWork = std::min(
      absoluteMaxRefinementWork, result.values.size() * std::size_t{64});
  if (boost::num_edges(g1) > maxCachedEdges ||
      boost::num_edges(g2) > maxCachedEdges) {
    return result;
  }
  const auto queryIncidentEdges = BuildIncidentEdgeCache(g1);
  const auto targetIncidentEdges = BuildIncidentEdgeCache(g2);
  NeighborhoodWorkspace workspace;
  std::size_t work = 0;
  for (unsigned int pass = 0; pass < maxRefinementPasses; ++pass) {
    bool changed = false;
    for (node_id query = 0; query < querySize; ++query) {
      for (node_id target = 0; target < targetSize; ++target) {
        if (!result.get(query, target)) {
          continue;
        }
        const auto matchResult = HasCompatibleNeighborhood(
            queryIncidentEdges[query], targetIncidentEdges[target], result,
            edgeCompatible, work, maxRefinementWork, workspace);
        if (matchResult == NeighborhoodMatchResult::Incomplete ||
            RDKit::ControlCHandler::getGotSignal()) {
          return result;
        }
        if (matchResult == NeighborhoodMatchResult::NoMatch) {
          result.remove(query, target);
          changed = true;
          if (!result.rowCounts[query]) {
            return result;
          }
        }
      }
    }
    if (!changed) {
      break;
    }
  }
  return result;
}

/*----------------------------------------------------
 * Sorts the nodes of a graphs, returning a
 * heap-allocated vector (using new) with the node ids
 * in the proper orders.
 * The sorting criterion takes into account:
 *    1 - The number of nodes with the same in/out
 *        degree.
 *    2 - The valence of the nodes.
 *    3 - The number of compatible target candidates, when available.
 *    4 - The number of already ordered neighbors.
 * The nodes at the beginning of the vector are
 * the most singular, from which the matching should start. Components are
 * kept contiguous and ring closures are prioritized.
 *--------------------------------------------------*/
template <class Graph>
node_id *SortNodesByFrequency(
    const Graph *g, const CandidateMatrix *candidateMatrix = nullptr) {
  std::vector<NodeInfo> vect;
  vect.reserve(boost::num_vertices(*g));
  typename Graph::vertex_iterator bNode, eNode;
  boost::tie(bNode, eNode) = boost::vertices(*g);
  while (bNode != eNode) {
    NodeInfo t;
    t.id = vect.size();
    t.in = t.out =
        boost::out_degree(*bNode, *g);  // <- assuming undirected graph
    t.candidates = candidateMatrix ? candidateMatrix->rowCounts[t.id] : 0;
    vect.push_back(t);
    ++bNode;
  }
  std::sort(vect.begin(), vect.end(), nodeInfoComp1);

  unsigned int run = 1;
  for (unsigned int i = 0; i < vect.size(); i += run) {
    for (run = 1; i + run < vect.size() && vect[i + run].in == vect[i].in &&
                  vect[i + run].out == vect[i].out;
         ++run) {
      ;
    }
    for (unsigned int j = 0; j < run; ++j) {
      vect[i + j].in += vect[i + j].out;
      vect[i + j].out = run;
    }
  }
  std::sort(vect.begin(), vect.end(), nodeInfoComp2);

  struct Candidate {
    node_id id;
    unsigned int selectedNbrs;
    unsigned int rank;
  };
  auto candidateComp = [](const Candidate &a, const Candidate &b) {
    if (a.selectedNbrs != b.selectedNbrs) {
      return a.selectedNbrs < b.selectedNbrs;
    }
    return a.rank > b.rank;
  };

  // A linear scan has lower overhead for normal-sized queries, while the
  // lazy heap prevents ordering itself from becoming quadratic for very
  // large queries. Both paths produce the same deterministic order.
  constexpr unsigned int heapThreshold = 128;
  const bool useHeap = vect.size() >= heapThreshold;
  std::priority_queue<Candidate, std::vector<Candidate>,
                      decltype(candidateComp)>
      candidates(candidateComp);
  std::vector<unsigned int> ranks;
  std::vector<unsigned char> selected(vect.size(), false);
  std::vector<unsigned int> selectedNbrs(vect.size(), 0);
  auto nodes = std::make_unique<node_id[]>(vect.size());

  if (useHeap) {
    ranks.resize(vect.size());
    for (unsigned int rank = 0; rank < vect.size(); ++rank) {
      const auto id = vect[rank].id;
      ranks[id] = rank;
      candidates.push({id, 0, rank});
    }
  }

  for (unsigned int i = 0; i < vect.size(); ++i) {
    node_id best = NULL_NODE;
    if (useHeap) {
      while (selected[candidates.top().id] ||
             candidates.top().selectedNbrs !=
                 selectedNbrs[candidates.top().id]) {
        candidates.pop();
      }
      best = candidates.top().id;
      candidates.pop();
    } else {
      unsigned int mostSelectedNbrs = 0;
      for (const auto &node : vect) {
        if (!selected[node.id] &&
            (best == NULL_NODE || selectedNbrs[node.id] > mostSelectedNbrs)) {
          best = node.id;
          mostSelectedNbrs = selectedNbrs[node.id];
        }
      }
    }
    assert(best != NULL_NODE);
    nodes[i] = best;
    selected[best] = true;

    typename Graph::adjacency_iterator nbrBeg, nbrEnd;
    boost::tie(nbrBeg, nbrEnd) = boost::adjacent_vertices(best, *g);
    while (nbrBeg != nbrEnd) {
      if (!selected[*nbrBeg]) {
        ++selectedNbrs[*nbrBeg];
        if (useHeap) {
          candidates.push({static_cast<node_id>(*nbrBeg), selectedNbrs[*nbrBeg],
                           ranks[*nbrBeg]});
        }
      }
      ++nbrBeg;
    }
  }

  return nodes.release();
}

/*----------------------------------------------------------
 * class VF2SubState
 * A representation of the SSS current state
 ---------------------------------------------------------*/
template <class Graph, class VertexCompatible, class EdgeCompatible,
          class MatchChecking>
class VF2SubState {
 private:
  Graph *g1, *g2;
  VertexCompatible &vc;
  EdgeCompatible &ec;
  MatchChecking &mc;
  unsigned int n1, n2;
  const CandidateMatrix *candidateMatrix;
  bool exactMatch;

  unsigned int core_len;
  unsigned int t1_len;
  unsigned int t2_len;  // Core nodes are also counted by these...
  node_id *core_1;
  node_id *core_2;
  node_id *term_1;
  node_id *term_2;

  node_id *order;

  long *share_count;
  int *vs_compared;

 public:
  VF2SubState(Graph *ag1, Graph *ag2, VertexCompatible &avc,
              EdgeCompatible &aec, MatchChecking &amc, bool sortNodes = false,
              const CandidateMatrix *candidates = nullptr,
              bool requireExactMatch = false)
      : g1(ag1),
        g2(ag2),
        vc(avc),
        ec(aec),
        mc(amc),
        n1(num_vertices(*ag1)),
        n2(num_vertices(*ag2)),
        candidateMatrix(candidates),
        exactMatch(requireExactMatch) {
    std::unique_ptr<node_id[]> newOrder;
    if (sortNodes) {
      newOrder.reset(SortNodesByFrequency(ag1, candidateMatrix));
    }

    core_len = 0;
    t1_len = 0;
    t2_len = 0;

    auto stateData =
        std::make_unique<node_id[]>(2 * (static_cast<std::size_t>(n1) + n2));
    auto newShareCount = std::make_unique<long>(1);
    core_1 = stateData.get();
    term_1 = core_1 + n1;
    core_2 = term_1 + n1;
    term_2 = core_2 + n2;

    for (unsigned int i = 0; i < n1; i++) {
      core_1[i] = NULL_NODE;
      term_1[i] = 0;
    }
    for (unsigned int i = 0; i < n2; i++) {
      core_2[i] = NULL_NODE;
      term_2[i] = 0;
    }
    vs_compared = nullptr;
    // vs_compared = new int[n1*n2];
    // memset((void *)vs_compared,0,n1*n2*sizeof(int));

    // es_compared = new std::map<unsigned int,bool>();
    order = newOrder.release();
    share_count = newShareCount.release();
    stateData.release();
  }

  VF2SubState(const VF2SubState &state)
      : g1(state.g1),
        g2(state.g2),
        vc(state.vc),
        ec(state.ec),
        mc(state.mc),
        n1(state.n1),
        n2(state.n2),
        candidateMatrix(state.candidateMatrix),
        exactMatch(state.exactMatch),
        order(state.order),
        vs_compared(state.vs_compared)
  // es_compared(state.es_compared)
  {
    core_len = state.core_len;
    t1_len = state.t1_len;
    t2_len = state.t2_len;

    core_1 = state.core_1;
    core_2 = state.core_2;
    term_1 = state.term_1;
    term_2 = state.term_2;
    share_count = state.share_count;

    ++(*share_count);
  }

  ~VF2SubState() {
    if (--*share_count == 0) {
      delete[] core_1;
      delete share_count;
      delete[] order;
      // delete [] vs_compared;
      // delete es_compared;
    }
  }

  bool IsGoal() { return core_len == n1; }
  bool MatchChecks(const node_id c1[], const node_id c2[]) {
    return mc(c1, c2);
  }
  bool IsDead() {
    return n1 > n2 || (exactMatch ? t1_len != t2_len : t1_len > t2_len);
  }
  unsigned int CoreLen() { return core_len; }
  Graph *GetGraph1() { return g1; }
  Graph *GetGraph2() { return g2; }

  bool NextPair(Pair<Graph> &pair) {
    if (pair.n1 == NULL_NODE) {
      pair.n1 = order == nullptr ? 0 : order[core_len];
    }
    if (pair.n2 == NULL_NODE) {
      pair.n2 = 0;
    } else {
      pair.n2++;
    }

#if 0
    std::cerr<<" **** np: "<< prev_n1<<","<<prev_n2<<std::endl;
    std::cerr<<"in_1 ";
    for(unsigned int i=0;i<n1;++i){
      std::cerr<<"("<<in_1[i]<<","<<out_1[i]<<"), ";
    } 
    std::cerr<<std::endl;
    std::cerr<<"in_2 ";
    for(unsigned int i=0;i<n2;++i){
      std::cerr<<"("<<in_2[i]<<","<<out_2[i]<<"), ";
    } 
    std::cerr<<std::endl;
#endif
    if (t1_len > core_len && t2_len > core_len) {
      if (order == nullptr) {
        while (pair.n1 < n1 &&
               (core_1[pair.n1] != NULL_NODE || term_1[pair.n1] == 0)) {
          pair.n1++;
          pair.n2 = 0;
        }
      } else {
        assert(core_1[pair.n1] == NULL_NODE && term_1[pair.n1] != 0);
      }

      /* Initialize VF2 Plus neighbor iterator.
       * The next query node (pair.n1) has been selected from the terminal
       * set and is therefore adjacent to an already mapped atom (in
       * core_1). Rather than select pair.n2 from all atoms (0...n2) we can
       * select it from the neighbors of this mapped atom (0...deg(nbor))
       * since it must also be adajcent to this mapped atom!
       */
      if (!pair.hasiter) {
        RDK_ADJ_ITER n1iter_beg, n1iter_end;
        boost::tie(n1iter_beg, n1iter_end) =
            boost::adjacent_vertices(pair.n1, *g1);

        node_id bestTargetNbr = NULL_NODE;
        unsigned int mappedNbrCount = 0;
        while (n1iter_beg != n1iter_end) {
          if (core_1[*n1iter_beg] != NULL_NODE) {
            if (bestTargetNbr == NULL_NODE) {
              bestTargetNbr = core_1[*n1iter_beg];
            }
            ++mappedNbrCount;
          }
          ++n1iter_beg;
        }

        assert(bestTargetNbr != NULL_NODE);

        // With multiple mapped query neighbors, optimized searches enumerate
        // candidates from the target neighborhood containing the fewest
        // currently unmapped atoms. Every valid candidate must occur in all of
        // those neighborhoods. Legacy searches retain their original neighbor
        // choice because changing it can change the order of returned matches.
        if (order != nullptr && mappedNbrCount > 1) {
          bestTargetNbr = NULL_NODE;
          unsigned int fewestCandidates = 0;
          boost::tie(n1iter_beg, n1iter_end) =
              boost::adjacent_vertices(pair.n1, *g1);
          while (n1iter_beg != n1iter_end) {
            if (core_1[*n1iter_beg] != NULL_NODE) {
              const auto targetNbr = core_1[*n1iter_beg];
              RDK_ADJ_ITER targetNbrBeg, targetNbrEnd;
              boost::tie(targetNbrBeg, targetNbrEnd) =
                  boost::adjacent_vertices(targetNbr, *g2);
              unsigned int candidateCount = 0;
              while (targetNbrBeg != targetNbrEnd &&
                     (bestTargetNbr == NULL_NODE ||
                      candidateCount < fewestCandidates)) {
                if (core_2[*targetNbrBeg] == NULL_NODE) {
                  ++candidateCount;
                }
                ++targetNbrBeg;
              }
              if (bestTargetNbr == NULL_NODE ||
                  candidateCount < fewestCandidates) {
                bestTargetNbr = targetNbr;
                fewestCandidates = candidateCount;
                if (!fewestCandidates) {
                  break;
                }
              }
            }
            ++n1iter_beg;
          }
        }

        boost::tie(pair.nbrbeg, pair.nbrend) =
            boost::adjacent_vertices(bestTargetNbr, *g2);
        pair.hasiter = true;
      }
    } else if (order == nullptr) {
      while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
        pair.n1++;
        pair.n2 = 0;
      }
    } else {
      assert(core_1[pair.n1] == NULL_NODE);
    }

    /* VF2 Plus iterator available? */
    if (pair.hasiter) {
      while (
          pair.nbrbeg < pair.nbrend &&
          (core_2[*pair.nbrbeg] != NULL_NODE ||
           (candidateMatrix && !candidateMatrix->get(pair.n1, *pair.nbrbeg)))) {
        ++pair.nbrbeg;
      }

      if (pair.nbrbeg < pair.nbrend) {
        pair.n2 = *pair.nbrbeg;
        ++pair.nbrbeg;
      } else {
        pair.n2 = n2;
      }
    } else if (t1_len > core_len && t2_len > core_len) {
      while (pair.n2 < n2 &&
             (core_2[pair.n2] != NULL_NODE || term_2[pair.n2] == 0 ||
              (candidateMatrix && !candidateMatrix->get(pair.n1, pair.n2)))) {
        pair.n2++;
      }
    } else {
      while (pair.n2 < n2 &&
             (core_2[pair.n2] != NULL_NODE ||
              (candidateMatrix && !candidateMatrix->get(pair.n1, pair.n2)))) {
        pair.n2++;
      }
    }
    return pair.n1 < n1 && pair.n2 < n2;
  }
  bool IsFeasiblePair(node_id node1, node_id node2) {
    assert(node1 < n1);
    assert(node2 < n2);
    assert(core_1[node1] == NULL_NODE);
    assert(core_2[node2] == NULL_NODE);

    // std::cerr<<"  ifp:"<<node1<<"-"<<node2<<"
    // "<<vs_compared->size()<<std::endl;
    // int &isCompat=vs_compared[node1*n2+node2];
    // if(isCompat==0){
    //   isCompat=vc(node1,node2)?1:-1;
    // }
    // if( isCompat<0 ){
    //   //std::cerr<<"  short1"<<std::endl;
    //   return false;
    // }

    const auto queryDegree = boost::out_degree(node1, *g1);
    const auto targetDegree = boost::out_degree(node2, *g2);
    if (candidateMatrix) {
      if (!candidateMatrix->get(node1, node2)) {
        return false;
      }
    } else {
      const bool degreeMatches = exactMatch ? queryDegree == targetDegree
                                            : queryDegree <= targetDegree;
      if (!degreeMatches || !vc(node1, node2)) {
        return false;
      }
    }

    unsigned int other1, other2;
    const bool useTerminalPruning = candidateMatrix || exactMatch;
    unsigned int term1 = 0, term2 = 0;
    unsigned int new1 = 0, new2 = 0;

    // Check the out edges of node1
    typename Graph::out_edge_iterator bNbrs, eNbrs;
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    while (bNbrs != eNbrs) {
      other1 = getOtherIdx(*g1, *bNbrs, node1);
      if (core_1[other1] != NULL_NODE) {
        other2 = core_1[other1];
        typename Graph::edge_descriptor oEdge;
        bool found;
        boost::tie(oEdge, found) = boost::edge(node2, other2, *g2);
        if (!found || !ec(*bNbrs, oEdge)) {
          // std::cerr<<"  short2"<<std::endl;
          return false;
        }
      } else if (useTerminalPruning) {
        if (term_1[other1]) ++term1;
        if (!term_1[other1]) ++new1;
      }
      ++bNbrs;
    }

    if (useTerminalPruning) {
      // Count target frontier neighbors. For an exact graph match, also reject
      // target edges to mapped nodes which have no corresponding query edge.
      boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
      while (bNbrs != eNbrs) {
        other2 = getOtherIdx(*g2, *bNbrs, node2);
        if (core_2[other2] != NULL_NODE) {
          if (exactMatch) {
            other1 = core_2[other2];
            // Compatibility for this edge was already checked while scanning
            // query edges above; this reverse check rejects extra target edges.
            if (!boost::edge(node1, other1, *g1).second) {
              return false;
            }
          }
        } else {
          if (term_2[other2]) ++term2;
          if (!term_2[other2]) ++new2;
        }
        ++bNbrs;
      }

      if (exactMatch) {
        return term1 == term2 && (term1 + new1) == (term2 + new2);
      }
      return term1 <= term2 && (term1 + new1) <= (term2 + new2);
    }
    return true;
  }
  void AddPair(node_id node1, node_id node2) {
    assert(node1 < n1);
    assert(node2 < n2);
    assert(core_len < n1);
    assert(core_len < n2);

    ++core_len;
    if (!term_1[node1]) {
      term_1[node1] = core_len;
      ++t1_len;
    }

    if (!term_2[node2]) {
      term_2[node2] = core_len;
      ++t2_len;
    }

    core_1[node1] = node2;
    core_2[node2] = node1;

    typename Graph::out_edge_iterator bNbrs, eNbrs;
    // FIX: this is explicitly ignoring directionality
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
      if (!term_1[other]) {
        term_1[other] = core_len;
        ++t1_len;
      }
      ++bNbrs;
    }

    // FIX: this is explicitly ignoring directionality
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
      if (!term_2[other]) {
        term_2[other] = core_len;
        ++t2_len;
      }
      ++bNbrs;
    }
  }
  void GetCoreSet(node_id c1[], node_id c2[]) {
    for (unsigned int i = 0; i < n1; ++i) {
      c1[i] = i;
      c2[i] = core_1[i];
    }
  }
  VF2SubState *Clone() { return new VF2SubState(*this); }
  void BackTrack(node_id node1, node_id node2) {
    if (term_1[node1] == core_len) {
      term_1[node1] = 0;
      --t1_len;
    }

    typename Graph::out_edge_iterator bNbrs, eNbrs;
    boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
      if (term_1[other] == core_len) {
        term_1[other] = 0;
        --t1_len;
      }
      ++bNbrs;
    }

    if (term_2[node2] == core_len) {
      term_2[node2] = 0;
      --t2_len;
    }

    boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
    while (bNbrs != eNbrs) {
      unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
      if (term_2[other] == core_len) {
        term_2[other] = 0;
        --t2_len;
      }
      ++bNbrs;
    }

    core_1[node1] = NULL_NODE;
    core_2[node2] = NULL_NODE;
    --core_len;
  }
  bool Match(node_id c1[], node_id c2[]) {
    if (IsGoal()) {
      GetCoreSet(c1, c2);
      if (MatchChecks(c1, c2)) {
        return true;
      }
      return false;
    }

    if (IsDead()) {
      return false;
    }

    Pair<Graph> pair;
    while (NextPair(pair) && !RDKit::ControlCHandler::getGotSignal()) {
      if (IsFeasiblePair(pair.n1, pair.n2)) {
        AddPair(pair.n1, pair.n2);
        if (Match(c1, c2)) {  // recurse
          return true;
        }
        BackTrack(pair.n1, pair.n2);
      }
    }
    return false;
  }

  template <class DoubleBackInsertionSequence>
  bool MatchAll(node_id c1[], node_id c2[], DoubleBackInsertionSequence &res,
                unsigned int lim = 0) {
    if (IsGoal()) {
      GetCoreSet(c1, c2);
      if (MatchChecks(c1, c2)) {
        typename DoubleBackInsertionSequence::value_type newSeq;
        newSeq.reserve(core_len);
        for (unsigned int i = 0; i < core_len; ++i) {
          newSeq.emplace_back(c1[i], c2[i]);
        }
        res.push_back(newSeq);
        return lim && res.size() >= lim;
      }
      return false;
    }

    if (IsDead()) {
      return false;
    }

    Pair<Graph> pair;
    while (NextPair(pair) && !RDKit::ControlCHandler::getGotSignal()) {
      if (IsFeasiblePair(pair.n1, pair.n2)) {
        AddPair(pair.n1, pair.n2);
        if (MatchAll(c1, c2, res, lim)) {  // recurse
          return true;
        }
        BackTrack(pair.n1, pair.n2);
      }
    }
    return false;
  }
};

/*-------------------------------------------------------------
 * static bool match(pn, c1, c2, s)
 * Finds a matching between two graphs, if it exists, starting
 * from state s.
 * Returns true a match has been found.
 * *pn is assigned the numbero of matched nodes, and
 * c1 and c2 will contain the ids of the corresponding nodes
 * in the two graphs.
 ------------------------------------------------------------*/
template <class SubState>
bool match(int *pn, node_id c1[], node_id c2[], SubState &s) {
  if (s.Match(c1, c2)) {
    // not needed, pn = num query atoms (n1)...
    *pn = s.CoreLen();
    return true;
  }
  return false;
}

/*-------------------------------------------------------------
 * static bool match(c1, c2, vis, usr_data, pcount)
 * Visits all the matchings between two graphs,  starting
 * from state s.
 * Returns true if the caller must stop the visit.
 * Stops when max_results is reached, set max_results to 0 to
 * keep going until there are no more matches
 *
 ------------------------------------------------------------*/
template <class SubState, class DoubleBackInsertionSequence>
bool match(node_id c1[], node_id c2[], SubState &s,
           DoubleBackInsertionSequence &res, unsigned int max_results) {
  s.MatchAll(c1, c2, res, max_results);
  return !res.empty();
}

struct AcceptAllFinalChecker {
  bool operator()(const node_id[], const node_id[]) const { return true; }
};

// Run the optimized, label-compatible traversal without invoking the final
// match checker, whose uniqueness state must not be changed by the precheck. A
// negative result proves that the full search cannot match; a positive result
// is followed by the legacy traversal to preserve ordering.
template <class Graph, class VertexLabeling, class EdgeLabeling>
bool hasPotentialMatch(const Graph &g1, const Graph &g2,
                       VertexLabeling &vertex_labeling,
                       EdgeLabeling &edge_labeling,
                       const CandidateMatrix *candidateMatrix,
                       bool exactMatch) {
  AcceptAllFinalChecker matchChecker;
  VF2SubState<const Graph, VertexLabeling, EdgeLabeling, AcceptAllFinalChecker>
      state(&g1, &g2, vertex_labeling, edge_labeling, matchChecker, true,
            candidateMatrix, exactMatch);
  auto ni1 = std::make_unique<node_id[]>(num_vertices(g1));
  auto ni2 = std::make_unique<node_id[]>(num_vertices(g2));
  int n = 0;
  return match(&n, ni1.get(), ni2.get(), state);
}
};  // end of namespace detail

template <
    class Graph, class VertexLabeling  // binary predicate
    ,
    class EdgeLabeling  // binary predicate
    ,
    class MatchChecking  // binary predicate
    ,
    class
    BackInsertionSequence  // contains
                           // std::pair<vertex_descriptor,vertex_descriptor>
    >
bool vf2(const Graph &g1, const Graph &g2, VertexLabeling &vertex_labeling,
         EdgeLabeling &edge_labeling, MatchChecking &match_checking,
         BackInsertionSequence &F) {
  F.clear();
  if (num_vertices(g1) > num_vertices(g2) || num_edges(g1) > num_edges(g2)) {
    return false;
  }

  RDKit::ControlCHandler hdlr;

  const bool exactMatch =
      num_vertices(g1) == num_vertices(g2) && num_edges(g1) == num_edges(g2);
  auto candidateMatrix = detail::BuildCandidateMatrix(
      g1, g2, vertex_labeling, edge_labeling, exactMatch);
  const auto *candidateMatrixPtr = candidateMatrix ? &candidateMatrix : nullptr;

  // Use the optimized query order to reject impossible matches without
  // changing the order in which successful matches are returned.
  if ((!candidateMatrixPtr || !candidateMatrix.hasEmptyRow()) &&
      detail::hasPotentialMatch(g1, g2, vertex_labeling, edge_labeling,
                                candidateMatrixPtr, exactMatch)) {
    detail::VF2SubState<const Graph, VertexLabeling, EdgeLabeling,
                        MatchChecking>
        s0(&g1, &g2, vertex_labeling, edge_labeling, match_checking, false,
           candidateMatrixPtr, exactMatch);
    auto ni1 = std::make_unique<detail::node_id[]>(num_vertices(g1));
    auto ni2 = std::make_unique<detail::node_id[]>(num_vertices(g2));
    int n = 0;
    if (match(&n, ni1.get(), ni2.get(), s0)) {
      auto sz = num_vertices(g1);
      F.reserve(sz);
      for (unsigned int i = 0; i < sz; ++i) {
        F.emplace_back(ni1[i], ni2[i]);
      }
    }
  }

  if (hdlr.getGotSignal()) {
    BOOST_LOG(rdWarningLog)
        << "Substructure search was interrupted, result may not include all matches"
        << std::endl;
  }

  return !F.empty();
};
template <class Graph, class VertexLabeling  // binary predicate
          ,
          class EdgeLabeling  // binary predicate
          ,
          class MatchChecking  // binary predicate
          ,
          class DoubleBackInsertionSequence  // contains a back insertion
                                             // sequence
          >
bool vf2_all(const Graph &g1, const Graph &g2, VertexLabeling &vertex_labeling,
             EdgeLabeling &edge_labeling, MatchChecking &match_checking,
             DoubleBackInsertionSequence &F, unsigned int max_results = 1000) {
  F.clear();
  if (num_vertices(g1) > num_vertices(g2) || num_edges(g1) > num_edges(g2)) {
    return false;
  }

  RDKit::ControlCHandler hdlr;

  const bool exactMatch =
      num_vertices(g1) == num_vertices(g2) && num_edges(g1) == num_edges(g2);
  auto candidateMatrix = detail::BuildCandidateMatrix(
      g1, g2, vertex_labeling, edge_labeling, exactMatch);
  const auto *candidateMatrixPtr = candidateMatrix ? &candidateMatrix : nullptr;

  // Use the optimized query order to reject impossible matches without
  // changing the order in which successful matches are returned.
  if ((!candidateMatrixPtr || !candidateMatrix.hasEmptyRow()) &&
      detail::hasPotentialMatch(g1, g2, vertex_labeling, edge_labeling,
                                candidateMatrixPtr, exactMatch)) {
    detail::VF2SubState<const Graph, VertexLabeling, EdgeLabeling,
                        MatchChecking>
        s0(&g1, &g2, vertex_labeling, edge_labeling, match_checking, false,
           candidateMatrixPtr, exactMatch);
    auto ni1 = std::make_unique<detail::node_id[]>(num_vertices(g1));
    auto ni2 = std::make_unique<detail::node_id[]>(num_vertices(g2));
    match(ni1.get(), ni2.get(), s0, F, max_results);
  }

  if (hdlr.getGotSignal()) {
    BOOST_LOG(rdWarningLog)
        << "Substructure search was interrupted, result may not include all matches"
        << std::endl;
  }

  return !F.empty();
};
}  // end of namespace boost
#endif

#undef RDK_ADJ_ITER
