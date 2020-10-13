#include "MUQ/Modeling/WorkGraphNode.h"

using namespace muq::Modeling;

WorkGraphNode::WorkGraphNode(std::shared_ptr<WorkPiece> piece, std::string const& name) : piece(piece), name(name) {
  assert(piece);
}
