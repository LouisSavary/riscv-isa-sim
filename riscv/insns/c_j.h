require_extension('C');
set_pc(pc + insn.rvc_j_imm());
p->update_predictor(pc, true, OpType::OPTYPE_JMP_DIRECT_UNCOND, npc);
