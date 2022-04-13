require_extension('C');
require(insn.rvc_rs1() != 0);
set_pc(RVC_RS1 & ~reg_t(1));

p->update_predictor(pc, true, OpType::OPTYPE_JMP_INDIRECT_UNCOND, npc,insn.bits());