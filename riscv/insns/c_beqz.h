require_extension('C');
bool taken = (RVC_RS1S == 0);
p->get_prediction(pc);
if (taken)
  set_pc(pc + insn.rvc_b_imm());
p->update_predictor(pc, taken, OpType::OPTYPE_JMP_DIRECT_COND, taken?npc:pc+4,insn.bits());