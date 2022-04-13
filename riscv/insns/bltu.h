bool taken = (RS1 < RS2);
p->get_prediction(pc);
if(taken)
  set_pc(BRANCH_TARGET);
p->update_predictor(pc, taken, OpType::OPTYPE_JMP_DIRECT_COND, taken?npc:pc+4,insn.bits());