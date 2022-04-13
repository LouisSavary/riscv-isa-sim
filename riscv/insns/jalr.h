reg_t tmp = npc;
set_pc((RS1 + insn.i_imm()) & ~reg_t(1));
WRITE_RD(tmp);

p->update_predictor(pc, true, OpType::OPTYPE_JMP_INDIRECT_UNCOND, npc,insn.bits());