require_extension('C');
if (xlen == 32) {
  reg_t tmp = npc;
  set_pc(pc + insn.rvc_j_imm());
  WRITE_REG(X_RA, tmp);
  p->update_predictor(pc, true, OpType::OPTYPE_JMP_DIRECT_UNCOND, npc);
} else { // c.addiw
  require(insn.rvc_rd() != 0);
  WRITE_RD(sext32(RVC_RS1 + insn.rvc_imm()));
}