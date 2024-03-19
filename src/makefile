localbin = ../bin

default: pic

pic: pic.o \
		common.o \
		moments.o \
		vector-math.o \
		file-io.o \
		lhs.o \
		logging.o \
		parameters.o \
		particles.o \
		positions.o \
		potential.o \
		random.o \
		rhs.o \
		setup.o \
		velocities.o
	-$(CLINKER) -o ${localbin}/$@ $^ $(LDLIBS)

solver: solver.o \
		common.o \
		moments.o \
		vector-math.o \
		file-io.o \
		lhs.o \
		logging.o \
		parameters.o \
		particles.o \
		positions.o \
		potential.o \
		random.o \
		rhs.o \
		setup.o \
		velocities.o
	-$(CLINKER) -o ${localbin}/$@ $^ $(LDLIBS) $(SLEPC_EPS_LIB)

include ${PETSC_DIR}/lib/petsc/conf/variables
# The following line will implicitly include PETSc rules. If we decide to remove
# it, we will need to explicitly include `${PETSC_DIR}/lib/petsc/conf/rules`.
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

