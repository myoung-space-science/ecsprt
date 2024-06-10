localbin = ../bin

default: pic

all: pic solver

pic: pic.o \
		boundaries.o \
		collisions.o \
		common.o \
		moments.o \
		vector-math.o \
		file-io.o \
		lhs.o \
		logging.o \
		parameters.o \
		positions.o \
		potential.o \
		random.o \
		rhs.o \
		setup.o \
		velocities.o
	-$(CLINKER) -o ${localbin}/$@ $^ $(LDLIBS)

solver: solver.o \
		boundaries.o \
		collisions.o \
		common.o \
		moments.o \
		vector-math.o \
		file-io.o \
		lhs.o \
		logging.o \
		parameters.o \
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
ifdef SLEPC_DIR
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
else
include ${PETSC_DIR}/lib/petsc/conf/rules
endif

