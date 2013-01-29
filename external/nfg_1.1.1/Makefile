BIN = bin
SRC = phantom
CC = gcc

GSL_PATH = /usr/local

CFLAGS = -c -Wall -O2 -I./  -I$(GSL_PATH)/include/
LDFLAGS = -L$(GSL_PATH)/lib
LIBRARIES = -lgsl -lgslcblas -lm

SHARED_INCLUDES = $(SRC)/shared/control_point.h $(SRC)/shared/segment.h $(SRC)/shared/shared.h $(SRC)/shared/strand.h $(SRC)/shared/strand_collection.h $(SRC)/shared/isotropic_region.h $(SRC)/shared/bundle.h

SHARED_OBJECTS = $(SRC)/shared/control_point.o $(SRC)/shared/segment.o $(SRC)/shared/shared.o $(SRC)/shared/strand.o $(SRC)/shared/strand_collection.o $(SRC)/shared/isotropic_region.o $(SRC)/shared/bundle.o


all: $(BIN)/rand_init $(BIN)/optimise $(BIN)/subdiv $(BIN)/resample $(BIN)/trim $(BIN)/mri_sim $(BIN)/noisify $(BIN)/draw_rois


$(BIN)/rand_init : $(SHARED_OBJECTS) $(SRC)/rand_init/main.o $(SRC)/rand_init/rand_init.o 
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/rand_init/main.o $(SRC)/rand_init/rand_init.o -o $@ $(LIBRARIES)

$(BIN)/optimise : $(SHARED_OBJECTS) $(SRC)/optimise/main.o $(SRC)/optimise/cost_function.o $(SRC)/optimise/reference_block.o $(SRC)/optimise/optimise.o $(SRC)/optimise/sample_block.o $(SRC)/optimise/sample.o
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/optimise/main.o $(SRC)/optimise/sample.o $(SRC)/optimise/cost_function.o $(SRC)/optimise/reference_block.o $(SRC)/optimise/optimise.o $(SRC)/optimise/sample_block.o -o $@ $(LIBRARIES)

$(BIN)/subdiv : $(SHARED_OBJECTS) $(SRC)/subdiv/main.o $(SRC)/subdiv/subdiv.o 
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/subdiv/main.o $(SRC)/subdiv/subdiv.o -o $@ $(LIBRARIES)

$(BIN)/resample : $(SHARED_OBJECTS) $(SRC)/resample/main.o $(SRC)/resample/resample.o 
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/resample/main.o $(SRC)/resample/resample.o -o $@ $(LIBRARIES)

$(BIN)/trim : $(SHARED_OBJECTS) $(SRC)/trim/main.o $(SRC)/trim/trim.o $(SRC)/trim/strand_section.o
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/trim/main.o $(SRC)/trim/trim.o $(SRC)/trim/strand_section.o -o $@ $(LIBRARIES)


$(BIN)/mri_sim : $(SHARED_OBJECTS) $(SRC)/mri_sim/main.o $(SRC)/mri_sim/mri_sim.o $(SRC)/mri_sim/segment_stats.o $(SRC)/mri_sim/subvoxel.o $(SRC)/mri_sim/overlap_strands.o $(SRC)/mri_sim/sim_voxel_intensities.o $(SRC)/mri_sim/voxel.o $(SRC)/mri_sim/segment_register.o $(SRC)/mri_sim/strand_collection_stats.o $(SRC)/mri_sim/isotropic_region_register.o
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/mri_sim/main.o $(SRC)/mri_sim/mri_sim.o $(SRC)/mri_sim/segment_stats.o $(SRC)/mri_sim/subvoxel.o $(SRC)/mri_sim/overlap_strands.o $(SRC)/mri_sim/sim_voxel_intensities.o $(SRC)/mri_sim/voxel.o $(SRC)/mri_sim/segment_register.o $(SRC)/mri_sim/strand_collection_stats.o $(SRC)/mri_sim/isotropic_region_register.o -o $@ $(LIBRARIES)
 
$(BIN)/noisify : $(SHARED_OBJECTS) $(SRC)/noisify/main.o $(SRC)/noisify/noisify.o
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/noisify/main.o $(SRC)/noisify/noisify.o -o $@ $(LIBRARIES)


$(BIN)/draw_rois : $(SHARED_OBJECTS) $(SRC)/draw_rois/main.o $(SRC)/draw_rois/draw_rois.o
	$(CC) $(LDFLAGS) $(SHARED_OBJECTS) $(SRC)/draw_rois/main.o $(SRC)/draw_rois/draw_rois.o -o $@ $(LIBRARIES)


.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm $(BIN)/*
	-rm $(SRC)/shared/*.o
	-rm $(SRC)/mri_sim/*.o
	-rm $(SRC)/rand_init/*.o
	-rm $(SRC)/trim/*.o
	-rm $(SRC)/subdiv/*.o
	-rm $(SRC)/resample/*.o
	-rm $(SRC)/optimise/*.o
	-rm $(SRC)/noisify/*.o
	-rm $(SRC)/draw_rois/*.o
