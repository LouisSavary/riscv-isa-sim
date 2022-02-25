reg_t tmp = npc;
set_pc(JUMP_TARGET);
WRITE_RD(tmp);

p->update_predictor(pc, true, OpType::OPTYPE_JMP_DIRECT_UNCOND, npc);