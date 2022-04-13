bool taken = (sreg_t(RS1) >= sreg_t(RS2));
p->get_prediction(pc);
if(taken)
  set_pc(BRANCH_TARGET);
p->update_predictor(pc, taken, OpType::OPTYPE_JMP_DIRECT_COND, npc,insn.bits());