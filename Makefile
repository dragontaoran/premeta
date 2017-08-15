CC = g++
EXE = premeta
CFLAGS = -O4 -std=c++11 -I./CodeBase -I./eigen-eigen-07105f7124f9
CODEBASE_OBJS = ./CodeBase/obj/command_line_utils.o ./CodeBase/obj/read_file_utils.o \
./CodeBase/obj/map_utils.o ./CodeBase/obj/data_structures.o \
./CodeBase/obj/eq_solver.o ./CodeBase/obj/cdf_fns.o ./CodeBase/obj/number_comparison.o ./CodeBase/obj/constants.o \
./CodeBase/obj/gamma_fns.o ./CodeBase/obj/string_utils.o ./CodeBase/obj/test_utils.o ./CodeBase/obj/csv_utils.o ./CodeBase/obj/vcf_utils.o
PREMETA_DIR = ./CodeBase/PreMeta/obj
PREMETA_OBJS = ${PREMETA_DIR}/premeta_constants.o ${PREMETA_DIR}/premeta_utils.o ${PREMETA_DIR}/rm_qc.o ${PREMETA_DIR}/write_mass_utils.o \
${PREMETA_DIR}/write_metaskat_utils.o ${PREMETA_DIR}/write_seqmeta_utils.o ${PREMETA_DIR}/write_r_binary_utils.o ${PREMETA_DIR}/write_raremetal_utils.o \
${PREMETA_DIR}/read_mass_utils.o ${PREMETA_DIR}/read_metaskat_utils.o ${PREMETA_DIR}/read_r_binary_utils.o ${PREMETA_DIR}/read_raremetal_utils.o \
${PREMETA_DIR}/read_seqmeta_utils.o ${PREMETA_DIR}/premeta_structures.o ${PREMETA_DIR}/snp_hwe.o ${PREMETA_DIR}/premeta.o

all:
	(cd ./CodeBase; make)
	(cd ./CodeBase/PreMeta; make)
	$(CC) $(CFlAGS) -o $(EXE) ${CODEBASE_OBJS} ${PREMETA_OBJS} -lm -lz

clean:
	(cd ./CodeBase; make clean)
	(cd ./CodeBase/PreMeta; make clean)
	rm -rf $(EXE)
