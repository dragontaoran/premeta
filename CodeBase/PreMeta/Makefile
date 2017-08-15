CC = g++
CFLAGS = -c -O4 -std=c++11 -I./ -I../../eigen-eigen-07105f7124f9 -I../
OBJ_DIR = obj
CodeBaseObj = ../obj

all: premeta premeta_utils rm_qc write_mass_utils write_metaskat_utils write_seqmeta_utils write_r_binary_utils write_raremetal_utils \
read_mass_utils read_metaskat_utils read_r_binary_utils read_raremetal_utils read_seqmeta_utils premeta_constants premeta_structures snp_hwe

premeta: ./premeta.cpp \
premeta_constants premeta_structures premeta_utils \
${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./premeta.cpp -o ./${OBJ_DIR}/premeta.o

premeta_utils: ./premeta_utils.cpp ./premeta_utils.h \
premeta_constants premeta_structures read_mass_utils read_metaskat_utils read_r_binary_utils read_raremetal_utils \
read_seqmeta_utils write_mass_utils write_metaskat_utils write_raremetal_utils write_r_binary_utils write_seqmeta_utils \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./premeta_utils.cpp -o ./${OBJ_DIR}/premeta_utils.o

rm_qc: ./rm_qc.cpp ./rm_qc.h \
premeta_constants premeta_structures read_raremetal_utils write_raremetal_utils \
${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./rm_qc.cpp -o ./${OBJ_DIR}/rm_qc.o

write_mass_utils: ./write_mass_utils.cpp ./write_mass_utils.h \
premeta_constants premeta_structures read_mass_utils \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./write_mass_utils.cpp -o ./${OBJ_DIR}/write_mass_utils.o

write_metaskat_utils: ./write_metaskat_utils.cpp ./write_metaskat_utils.h \
premeta_constants premeta_structures read_metaskat_utils \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./write_metaskat_utils.cpp -o ./${OBJ_DIR}/write_metaskat_utils.o

write_seqmeta_utils: ./write_seqmeta_utils.cpp ./write_seqmeta_utils.h \
premeta_constants premeta_structures read_raremetal_utils write_r_binary_utils \
${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./write_seqmeta_utils.cpp -o ./${OBJ_DIR}/write_seqmeta_utils.o

write_r_binary_utils: ./write_r_binary_utils.cpp ./write_r_binary_utils.h \
premeta_constants premeta_structures
	${CC} ${CFLAGS} ./write_r_binary_utils.cpp -o ./${OBJ_DIR}/write_r_binary_utils.o

write_raremetal_utils: ./write_raremetal_utils.cpp ./write_raremetal_utils.h \
premeta_constants premeta_structures snp_hwe \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./write_raremetal_utils.cpp -o ./${OBJ_DIR}/write_raremetal_utils.o

read_mass_utils: ./read_mass_utils.cpp ./read_mass_utils.h \
premeta_constants premeta_structures \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./read_mass_utils.cpp -o ./${OBJ_DIR}/read_mass_utils.o

read_metaskat_utils: ./read_metaskat_utils.cpp ./read_metaskat_utils.h \
premeta_constants premeta_structures \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./read_metaskat_utils.cpp -o ./${OBJ_DIR}/read_metaskat_utils.o

read_r_binary_utils: ./read_r_binary_utils.cpp ./read_r_binary_utils.h \
premeta_constants premeta_structures ${CodeBaseObj}/csv_utils.o \
${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o
	${CC} ${CFLAGS} ./read_r_binary_utils.cpp -o ./${OBJ_DIR}/read_r_binary_utils.o

read_raremetal_utils: ./read_raremetal_utils.cpp ./read_raremetal_utils.h \
premeta_constants premeta_structures \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./read_raremetal_utils.cpp -o ./${OBJ_DIR}/read_raremetal_utils.o

read_seqmeta_utils: ./read_seqmeta_utils.cpp ./read_seqmeta_utils.h \
premeta_constants premeta_structures \
${CodeBaseObj}/csv_utils.o ${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/number_comparison.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./read_seqmeta_utils.cpp -o ./${OBJ_DIR}/read_seqmeta_utils.o

premeta_constants: ./premeta_constants.cpp ./premeta_constants.h
	${CC} ${CFLAGS} ./premeta_constants.cpp -o ./${OBJ_DIR}/premeta_constants.o

premeta_structures: ./premeta_structures.cpp ./premeta_structures.h \
${CodeBaseObj}/vcf_utils.o ${CodeBaseObj}/map_utils.o ${CodeBaseObj}/string_utils.o
	${CC} ${CFLAGS} ./premeta_structures.cpp -o ./${OBJ_DIR}/premeta_structures.o

snp_hwe: ./snp_hwe.cpp ./snp_hwe.h
	${CC} ${CFLAGS} ./snp_hwe.cpp -o ./${OBJ_DIR}/snp_hwe.o

clean:
	rm -rf ./${OBJ_DIR}/*.o

