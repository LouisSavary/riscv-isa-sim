// See LICENSE for license details.
#ifndef _RISCV_PROCESSOR_H
#define _RISCV_PROCESSOR_H

#include "decode.h"
#include "config.h"
#include "trap.h"
#include "abstract_device.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <cassert>
#include "debug_rom_defines.h"
#include "entropy_source.h"
#include "csrs.h"
// #include <predictor.h>

#define NB_PRE_PRED 0

typedef enum {
  OPTYPE_OP               =2,

  OPTYPE_RET_UNCOND,
  OPTYPE_JMP_DIRECT_UNCOND,
  OPTYPE_JMP_INDIRECT_UNCOND,
  OPTYPE_CALL_DIRECT_UNCOND,
  OPTYPE_CALL_INDIRECT_UNCOND,

  OPTYPE_RET_COND,
  OPTYPE_JMP_DIRECT_COND,
  OPTYPE_JMP_INDIRECT_COND,
  OPTYPE_CALL_DIRECT_COND,
  OPTYPE_CALL_INDIRECT_COND,

  OPTYPE_ERROR,

  OPTYPE_MAX
}OpType;

class processor_t;
class mmu_t;
typedef reg_t (*insn_func_t)(processor_t*, insn_t, reg_t);
class simif_t;
class trap_t;
class extension_t;
class disassembler_t;
class control_graph;

reg_t illegal_instruction(processor_t* p, insn_t insn, reg_t pc);

struct insn_desc_t
{
  insn_bits_t match;
  insn_bits_t mask;
  insn_func_t rv32i;
  insn_func_t rv64i;
  insn_func_t rv32e;
  insn_func_t rv64e;

  insn_func_t func(int xlen, bool rve)
  {
    if (rve)
      return xlen == 64 ? rv64e : rv32e;
    else
      return xlen == 64 ? rv64i : rv32i;
  }

  static insn_desc_t illegal()
  {
    return {0, 0, &illegal_instruction, &illegal_instruction, &illegal_instruction, &illegal_instruction};
  }
};

// regnum, data
typedef std::unordered_map<reg_t, freg_t> commit_log_reg_t;

// addr, value, size
typedef std::vector<std::tuple<reg_t, uint64_t, uint8_t>> commit_log_mem_t;

typedef enum
{
  ACTION_DEBUG_EXCEPTION = MCONTROL_ACTION_DEBUG_EXCEPTION,
  ACTION_DEBUG_MODE = MCONTROL_ACTION_DEBUG_MODE,
  ACTION_TRACE_START = MCONTROL_ACTION_TRACE_START,
  ACTION_TRACE_STOP = MCONTROL_ACTION_TRACE_STOP,
  ACTION_TRACE_EMIT = MCONTROL_ACTION_TRACE_EMIT
} mcontrol_action_t;

typedef enum
{
  MATCH_EQUAL = MCONTROL_MATCH_EQUAL,
  MATCH_NAPOT = MCONTROL_MATCH_NAPOT,
  MATCH_GE = MCONTROL_MATCH_GE,
  MATCH_LT = MCONTROL_MATCH_LT,
  MATCH_MASK_LOW = MCONTROL_MATCH_MASK_LOW,
  MATCH_MASK_HIGH = MCONTROL_MATCH_MASK_HIGH
} mcontrol_match_t;

typedef struct
{
  uint8_t type;
  bool dmode;
  uint8_t maskmax;
  bool select;
  bool timing;
  mcontrol_action_t action;
  bool chain;
  mcontrol_match_t match;
  bool m;
  bool h;
  bool s;
  bool u;
  bool execute;
  bool store;
  bool load;
} mcontrol_t;

enum VRM{
  RNU = 0,
  RNE,
  RDN,
  ROD,
  INVALID_RM
};

template<uint64_t N>
struct type_usew_t;

template<>
struct type_usew_t<8>
{
  using type=uint8_t;
};

template<>
struct type_usew_t<16>
{
  using type=uint16_t;
};

template<>
struct type_usew_t<32>
{
  using type=uint32_t;
};

template<>
struct type_usew_t<64>
{
  using type=uint64_t;
};

template<uint64_t N>
struct type_sew_t;

template<>
struct type_sew_t<8>
{
  using type=int8_t;
};

template<>
struct type_sew_t<16>
{
  using type=int16_t;
};

template<>
struct type_sew_t<32>
{
  using type=int32_t;
};

template<>
struct type_sew_t<64>
{
  using type=int64_t;
};


// architectural state of a RISC-V hart
struct state_t
{
  void reset(processor_t* const proc, reg_t max_isa);

  static const int num_triggers = 4;

  reg_t pc;
  regfile_t<reg_t, NXPR, true> XPR;
  regfile_t<freg_t, NFPR, false> FPR;

  // control and status registers
  std::unordered_map<reg_t, csr_t_p> csrmap;
  reg_t prv;    // TODO: Can this be an enum instead?
  bool v;
  misa_csr_t_p misa;
  mstatus_csr_t_p mstatus;
  csr_t_p mepc;
  csr_t_p mtval;
  csr_t_p mtvec;
  csr_t_p mcause;
  minstret_csr_t_p minstret;
  mie_csr_t_p mie;
  mip_csr_t_p mip;
  csr_t_p medeleg;
  csr_t_p mideleg;
  csr_t_p mcounteren;
  csr_t_p scounteren;
  csr_t_p sepc;
  csr_t_p stval;
  csr_t_p stvec;
  virtualized_csr_t_p satp;
  csr_t_p scause;

  csr_t_p mtval2;
  csr_t_p mtinst;
  csr_t_p hstatus;
  csr_t_p hideleg;
  csr_t_p hedeleg;
  csr_t_p hcounteren;
  csr_t_p htval;
  csr_t_p htinst;
  csr_t_p hgatp;
  sstatus_csr_t_p sstatus;
  vsstatus_csr_t_p vsstatus;
  csr_t_p vstvec;
  csr_t_p vsepc;
  csr_t_p vscause;
  csr_t_p vstval;
  csr_t_p vsatp;

  csr_t_p dpc;
  dcsr_csr_t_p dcsr;
  csr_t_p tselect;
  mcontrol_t mcontrol[num_triggers];
  tdata2_csr_t_p tdata2;
  bool debug_mode;

  static const int max_pmp = 16;
  pmpaddr_csr_t_p pmpaddr[max_pmp];

  csr_t_p fflags;
  csr_t_p frm;

  csr_t_p menvcfg;
  csr_t_p senvcfg;
  csr_t_p henvcfg;

  bool serialized; // whether timer CSRs are in a well-defined state

  // When true, execute a single instruction and then enter debug mode.  This
  // can only be set by executing dret.
  enum {
      STEP_NONE,
      STEP_STEPPING,
      STEP_STEPPED
  } single_step;

#ifdef RISCV_ENABLE_COMMITLOG
  commit_log_reg_t log_reg_write;
  commit_log_mem_t log_mem_read;
  commit_log_mem_t log_mem_write;
  reg_t last_inst_priv;
  int last_inst_xlen;
  int last_inst_flen;
#endif
};

typedef enum {
  OPERATION_EXECUTE,
  OPERATION_STORE,
  OPERATION_LOAD,
} trigger_operation_t;

typedef enum {
  // 65('A') ~ 90('Z') is reserved for standard isa in misa
  EXT_ZFH,
  EXT_ZFHMIN,
  EXT_ZBA,
  EXT_ZBB,
  EXT_ZBC,
  EXT_ZBS,
  EXT_ZBKB,
  EXT_ZBKC,
  EXT_ZBKX,
  EXT_ZKND,
  EXT_ZKNE,
  EXT_ZKNH,
  EXT_ZKSED,
  EXT_ZKSH,
  EXT_ZKR,
  EXT_ZMMUL,
  EXT_ZBPBO,
  EXT_ZPN,
  EXT_ZPSFOPERAND,
  EXT_SVNAPOT,
  EXT_SVPBMT,
  EXT_SVINVAL,
  EXT_ZDINX,
  EXT_ZFINX,
  EXT_ZHINX,
  EXT_ZHINXMIN,
  EXT_XBITMANIP,
  EXT_ZICBOM,
  EXT_ZICBOZ,
} isa_extension_t;

typedef enum {
  IMPL_MMU_SV32,
  IMPL_MMU_SV39,
  IMPL_MMU_SV48,
  IMPL_MMU_SBARE,
  IMPL_MMU,
} impl_extension_t;

// Count number of contiguous 1 bits starting from the LSB.
static int cto(reg_t val)
{
  int res = 0;
  while ((val & 1) == 1)
    val >>= 1, res++;
  return res;
}

class isa_parser_t {
public:
  isa_parser_t(const char* str);
  ~isa_parser_t(){};
  unsigned get_max_xlen() { return max_xlen; }
  std::string get_isa_string() { return isa_string; }
  bool extension_enabled(unsigned char ext) const {
    if (ext >= 'A' && ext <= 'Z')
      return (max_isa >> (ext - 'A')) & 1;
    else
      return extension_table[ext];
  }
protected:
  unsigned max_xlen;
  reg_t max_isa;
  std::vector<bool> extension_table;
  std::string isa_string;
  std::unordered_map<std::string, extension_t*> custom_extensions;
};


typedef struct node_descr{
  uint32_t id = 0;
  OpType type = OPTYPE_JMP_INDIRECT_UNCOND;
  uint32_t opcode = 0;
  uint32_t taken_cnt = 0;
  uint32_t not_taken_cnt = 0;
  uint64_t pc = 0;
  uint32_t target_count = 0;
} node_t;

typedef struct edge_descr{
  uint32_t id=0;
  uint32_t size = 0;
  uint64_t target;
  uint64_t pc_src;
  uint64_t pc_dst;
  uint64_t taken_cnt;
  bool taken;
} edge_t;


class control_graph {
public: 

  control_graph():
      nodes(),
      edges(),
      edge_sequence() {
    branch(0, 0, OPTYPE_JMP_DIRECT_UNCOND, true, 0);
  }
  
  void dump(){

    FILE* f = stderr;//fopen("bt9_dump", "w");
    // if (f==nullptr) return;

    branch(-1, -1, OPTYPE_JMP_DIRECT_UNCOND, true, 0);

    //graph linearization, to remove interruptions hazards

    // display 
    uint64_t nb_nbr_instr=0;
    uint64_t nb_br_instr =0;

    auto it = edges.begin();
    while (it != edges.end()) {
      nb_nbr_instr += it->second.size;

      it++;
    }

    auto it2 = nodes.begin();
    while (it2 != nodes.end()) {
      nb_nbr_instr += it2->second.taken_cnt;
      nb_nbr_instr += it2->second.not_taken_cnt;

      it2++;
    }

    fprintf(f,"BT9_SPA_TRACE_FORMAT\n");
    fprintf(f,"bt9_minor_version: 0\n");
    fprintf(f,"has_physical_address: 0\n");
    fprintf(f,"md5_checksum:\n");
    fprintf(f,"conversion_date:\n");
    fprintf(f,"original_stf_input_file:\n");
    fprintf(f,"total_instruction_count: %16lu    # Instruction count\n", nb_nbr_instr+nb_br_instr);
    fprintf(f,"branch_instruction_count: %16lu    # Branch count\n", nb_br_instr);
    fprintf(f,"invalid_physical_branch_target_count:                0    # Invalid Physical Target Count \n");
    fprintf(f,"A32_instruction_count: %16lu    # A32 instructions \n", nb_nbr_instr+nb_br_instr);
    fprintf(f,"A64_instruction_count:                0    # A64 instructions\n");
    fprintf(f,"T32_instruction_count:                0    # T32 instructions\n");
    fprintf(f,"unidentified_instruction_count:                0    # Unidentified instructions\n");
    
    fprintf(f,"BT9_NODES\n#NODE    id      virtual_address    physical_address          opcode  size\n");


    auto it_node = nodes.begin();
    while (it_node != nodes.end()) {
      node_t n = it_node->second;
      if (n.id == id_node-1) {
        fprintf(f,"NODE %6u   0xffffffffffffffff                    -                0    0\n", id_node-1);
      } else if (n.id == 0) {
        fprintf(f,"NODE %6u   0x0000000000000000                    -                0    0\n", id_node-1);
      } else {
        const char* type = "";
        switch (n.type) {
        case OPTYPE_JMP_DIRECT_COND: type="JMP+DIR+CND"; break;
        case OPTYPE_JMP_DIRECT_UNCOND: type="JMP+DIR+UCD"; break;
        case OPTYPE_JMP_INDIRECT_COND: type="JMP+IND+CND"; break;
        case OPTYPE_JMP_INDIRECT_UNCOND: type="JMP+IND+UCD"; break;
        case OPTYPE_RET_UNCOND: type="RET+IND+UCD"; break;
        default: type="ERROR";
        }
        fprintf(f,"NODE %6u   0x%016lx                -       0x%08x   4  class:  %12s  behavior: DYN+DIR  taken_cnt:  %6u  not_taken_cnt:   %6u  tgt_cnt:   %u\n", 
            n.id, n.pc, n.opcode, type, n.taken_cnt, n.not_taken_cnt, n.target_count);
      }
      it_node++;
    }

    fprintf(f,"BT9_EDGES\n#EDGE    id  src_id   dest_id  taken      br_virt_target       br_phy_target      inst_cnt \n");
    auto it_edge = edges.begin();
    while (it_edge != edges.end()) {
      edge_t e = it_edge->second;
      fprintf(f,"EDGE  %6u  %6u  %6u       %s  0x%016lx             -        %10u    traverse_cnt:  %7lu\n", 
           e.id, nodes[e.pc_src].id, nodes[e.pc_dst].id, e.taken?"T":"N",e.target, e.size, e.taken_cnt);

      it_edge++;
    }

    fprintf(f,"BT9_EDGE_SEQUENCE\n");
    auto it_seq = edge_sequence.begin();
    while (it_seq != edge_sequence.end()) {
      fprintf(f,"%u\n", *(it_seq++));
    }
    // fclose(f);
  }


  void branch(uint64_t src, uint64_t dst, OpType op, bool taken, uint32_t opcode) {
    // printf("branch %16lx %16lx\n", src, dst);
    if (nodes.find(src) != nodes.end()) {
    //    if exists: increments
      nodes[src].taken_cnt += taken;
      nodes[src].not_taken_cnt += !taken;
    } 
    else {
    //    else: initialize it
      nodes[src].id = id_node++;
      nodes[src].pc = src;
      nodes[src].type = op;
      nodes[src].opcode = opcode;
      nodes[src].taken_cnt = taken;
      nodes[src].not_taken_cnt = !taken;
      nodes[src].target_count = 0;

      // order_nodes.push_back(nodes[src]);
    }


    // close current edge
    current_edge.pc_dst = src;
    uint64_t key = (long)((long)nodes[current_edge.pc_src].id)<<32 | (long)nodes[src].id;
    if (edges.find(key) != edges.end()){
    //    if exists: increments
      edges[key].taken_cnt ++;
    }
    else {
    //    else: create it
      // printf("0x%016lx\n", key);
      current_edge.taken_cnt = 1;
      current_edge.size = abs((long long signed)current_edge.pc_dst-(long long signed)current_edge.target)>>5;


      edges[key] = current_edge;

      nodes[current_edge.pc_src].target_count ++;
      // order_nodes.at(nodes[src].id).target_count++;


      current_edge.id = id_edge++;
    }
    edge_sequence.push_back(edges[key].id);

    // set new edge 
    current_edge.pc_src = src;
    current_edge.target = dst;
    current_edge.taken = taken;
    
  }

protected:
  std::unordered_map<uint64_t, node_t> nodes;
  // std::vector<node_t> order_nodes;
  std::unordered_map<uint64_t, edge_t> edges;
  std::vector<uint32_t> edge_sequence;
  edge_t current_edge;
  uint32_t id_node=0;
  uint32_t id_edge=0;
};

// this class represents one processor in a RISC-V machine.
class processor_t : public abstract_device_t, public isa_parser_t
{
public:
  processor_t(const char* isa, const char* priv, const char* varch,
              simif_t* sim, uint32_t id, bool halt_on_reset,
              FILE *log_file, std::ostream& sout_); // because of command line option --log and -s we need both
  ~processor_t();

  void set_debug(bool value);
  void set_histogram(bool value);
#ifdef RISCV_ENABLE_COMMITLOG
  void enable_log_commits();
  bool get_log_commits_enabled() const { return log_commits_enabled; }
#endif
  void reset();
  void step(size_t n); // run for n cycles
  void set_csr(int which, reg_t val);
  uint32_t get_id() const { return id; }
  reg_t get_csr(int which, insn_t insn, bool write, bool peek = 0);
  reg_t get_csr(int which) { return get_csr(which, insn_t(0), false, true); }
  mmu_t* get_mmu() { return mmu; }
  state_t* get_state() { return &state; }
  unsigned get_xlen() { return xlen; }
  unsigned get_const_xlen() {
    // Any code that assumes a const xlen should use this method to
    // document that assumption. If Spike ever changes to allow
    // variable xlen, this method should be removed.
    return xlen;
  }
  unsigned get_flen() {
    return extension_enabled('Q') ? 128 :
           extension_enabled('D') ? 64 :
           extension_enabled('F') ? 32 : 0;
  }
  extension_t* get_extension();
  extension_t* get_extension(const char* name);
  bool any_custom_extensions() const {
    return !custom_extensions.empty();
  }
  bool extension_enabled(unsigned char ext) const {
    if (ext >= 'A' && ext <= 'Z')
      return state.misa->extension_enabled(ext);
    else
      return extension_table[ext];
  }
  // Is this extension enabled? and abort if this extension can
  // possibly be disabled dynamically. Useful for documenting
  // assumptions about writable misa bits.
  bool extension_enabled_const(unsigned char ext) const {
    if (ext >= 'A' && ext <= 'Z')
      return state.misa->extension_enabled_const(ext);
    else
      return extension_table[ext];  // assume this can't change
  }
  void set_impl(uint8_t impl, bool val) { impl_table[impl] = val; }
  bool supports_impl(uint8_t impl) const {
    return impl_table[impl];
  }
  reg_t pc_alignment_mask() {
    return ~(reg_t)(extension_enabled('C') ? 0 : 2);
  }
  void check_pc_alignment(reg_t pc) {
    if (unlikely(pc & ~pc_alignment_mask()))
      throw trap_instruction_address_misaligned(state.v, pc, 0, 0);
  }
  reg_t legalize_privilege(reg_t);
  void set_privilege(reg_t);
  void set_virt(bool);
  void update_histogram(reg_t pc);
  const disassembler_t* get_disassembler() { return disassembler; }

  FILE *get_log_file() { return log_file; }

  void register_insn(insn_desc_t);
  void register_extension(extension_t*);

  // MMIO slave interface
  bool load(reg_t addr, size_t len, uint8_t* bytes);
  bool store(reg_t addr, size_t len, const uint8_t* bytes);

  // When true, display disassembly of each instruction that's executed.
  bool debug;
  // When true, take the slow simulation path.
  bool slow_path();
  bool halted() { return state.debug_mode; }
  enum {
    HR_NONE,    /* Halt request is inactive. */
    HR_REGULAR, /* Regular halt request/debug interrupt. */
    HR_GROUP    /* Halt requested due to halt group. */
  } halt_request;

  // Return the index of a trigger that matched, or -1.
  inline int trigger_match(trigger_operation_t operation, reg_t address, reg_t data)
  {
    if (state.debug_mode)
      return -1;

    bool chain_ok = true;

    for (unsigned int i = 0; i < state.num_triggers; i++) {
      if (!chain_ok) {
        chain_ok |= !state.mcontrol[i].chain;
        continue;
      }

      if ((operation == OPERATION_EXECUTE && !state.mcontrol[i].execute) ||
          (operation == OPERATION_STORE && !state.mcontrol[i].store) ||
          (operation == OPERATION_LOAD && !state.mcontrol[i].load) ||
          (state.prv == PRV_M && !state.mcontrol[i].m) ||
          (state.prv == PRV_S && !state.mcontrol[i].s) ||
          (state.prv == PRV_U && !state.mcontrol[i].u)) {
        continue;
      }

      reg_t value;
      if (state.mcontrol[i].select) {
        value = data;
      } else {
        value = address;
      }

      // We need this because in 32-bit mode sometimes the PC bits get sign
      // extended.
      if (xlen == 32) {
        value &= 0xffffffff;
      }

      auto tdata2 = state.tdata2->read(i);
      switch (state.mcontrol[i].match) {
        case MATCH_EQUAL:
          if (value != tdata2)
            continue;
          break;
        case MATCH_NAPOT:
          {
            reg_t mask = ~((1 << (cto(tdata2)+1)) - 1);
            if ((value & mask) != (tdata2 & mask))
              continue;
          }
          break;
        case MATCH_GE:
          if (value < tdata2)
            continue;
          break;
        case MATCH_LT:
          if (value >= tdata2)
            continue;
          break;
        case MATCH_MASK_LOW:
          {
            reg_t mask = tdata2 >> (xlen/2);
            if ((value & mask) != (tdata2 & mask))
              continue;
          }
          break;
        case MATCH_MASK_HIGH:
          {
            reg_t mask = tdata2 >> (xlen/2);
            if (((value >> (xlen/2)) & mask) != (tdata2 & mask))
              continue;
          }
          break;
      }

      if (!state.mcontrol[i].chain) {
        return i;
      }
      chain_ok = true;
    }
    return -1;
  }

  void trigger_updated();

  void set_pmp_num(reg_t pmp_num);
  void set_pmp_granularity(reg_t pmp_granularity);
  void set_mmu_capability(int cap);

  const char* get_symbol(uint64_t addr);

  bool get_prediction(uint64_t pc);
  // bool get_pre_prediction(uint64_t pc); sera fait dans get_pred
  void update_predictor(uint64_t pc, bool taken, OpType op, uint64_t target, uint32_t opcode);
  // uint* get_pred_stats();
  // uint64_t fetch_next_br(uint64_t pc);
  control_graph cg;

private:

  // OpType last_br_optype;
  // uint64_t last_br_target;

  // PREDICTOR* branch_predictor;
  // bool** predictions;
  // uint pred_id;
  // uint* pre_pred_stats;

  simif_t* sim;
  mmu_t* mmu; // main memory is always accessed via the mmu
  std::unordered_map<std::string, extension_t*> custom_extensions;
  disassembler_t* disassembler;
  state_t state;
  uint32_t id;
  unsigned xlen;
  bool histogram_enabled;
  bool log_commits_enabled;
  FILE *log_file;
  std::ostream sout_; // needed for socket command interface -s, also used for -d and -l, but not for --log
  bool halt_on_reset;
  std::vector<bool> impl_table;

  std::vector<insn_desc_t> instructions;
  std::map<reg_t,uint64_t> pc_histogram;

  static const size_t OPCODE_CACHE_SIZE = 8191;
  insn_desc_t opcode_cache[OPCODE_CACHE_SIZE];

  void take_pending_interrupt() { take_interrupt(state.mip->read() & state.mie->read()); }
  void take_interrupt(reg_t mask); // take first enabled interrupt in mask
  void take_trap(trap_t& t, reg_t epc); // take an exception
  void disasm(insn_t insn); // disassemble and print an instruction
  int paddr_bits();

  void enter_debug_mode(uint8_t cause);

  void debug_output_log(std::stringstream *s); // either output to interactive user or write to log file

  friend class mmu_t;
  friend class clint_t;
  friend class extension_t;

  void parse_varch_string(const char*);
  void parse_priv_string(const char*);
  void build_opcode_map();
  void register_base_instructions();
  insn_func_t decode_insn(insn_t insn);

  // Track repeated executions for processor_t::disasm()
  uint64_t last_pc, last_bits, executions;
public:
  entropy_source es; // Crypto ISE Entropy source.

  reg_t n_pmp;
  reg_t lg_pmp_granularity;
  reg_t pmp_tor_mask() { return -(reg_t(1) << (lg_pmp_granularity - PMP_SHIFT)); }

  class vectorUnit_t {
    public:
      processor_t* p;
      void *reg_file;
      char reg_referenced[NVPR];
      int setvl_count;
      reg_t vlmax;
      reg_t vlenb;
      csr_t_p vxsat;
      vector_csr_t_p vxrm, vstart, vl, vtype;
      reg_t vma, vta;
      reg_t vsew;
      float vflmul;
      reg_t ELEN, VLEN;
      bool vill;
      bool vstart_alu;

      // vector element for varies SEW
      template<class T>
        T& elt(reg_t vReg, reg_t n, bool is_write = false){
          assert(vsew != 0);
          assert((VLEN >> 3)/sizeof(T) > 0);
          reg_t elts_per_reg = (VLEN >> 3) / (sizeof(T));
          vReg += n / elts_per_reg;
          n = n % elts_per_reg;
#ifdef WORDS_BIGENDIAN
          // "V" spec 0.7.1 requires lower indices to map to lower significant
          // bits when changing SEW, thus we need to index from the end on BE.
          n ^= elts_per_reg - 1;
#endif
          reg_referenced[vReg] = 1;

#ifdef RISCV_ENABLE_COMMITLOG
          if (is_write)
            p->get_state()->log_reg_write[((vReg) << 4) | 2] = {0, 0};
#endif

          T *regStart = (T*)((char*)reg_file + vReg * (VLEN >> 3));
          return regStart[n];
        }
    public:

      void reset();

      vectorUnit_t():
        p(0),
        reg_file(0),
        reg_referenced{0},
        setvl_count(0),
        vlmax(0),
        vlenb(0),
        vxsat(0),
        vxrm(0),
        vstart(0),
        vl(0),
        vtype(0),
        vma(0),
        vta(0),
        vsew(0),
        vflmul(0),
        ELEN(0),
        VLEN(0),
        vill(false),
        vstart_alu(false) {
      }

      ~vectorUnit_t(){
        free(reg_file);
        reg_file = 0;
      }

      reg_t set_vl(int rd, int rs1, reg_t reqVL, reg_t newType);

      reg_t get_vlen() { return VLEN; }
      reg_t get_elen() { return ELEN; }
      reg_t get_slen() { return VLEN; }

      VRM get_vround_mode() {
        return (VRM)(vxrm->read());
      }
  };

  vectorUnit_t VU;
};



#endif
