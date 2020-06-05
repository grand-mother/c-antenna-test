CFLAGS += $(shell pkg-config --cflags hdf5) -O0 -g3
LIBS= $(shell pkg-config --libs hdf5)

# H5PY_LIBS= ../grand/bin/python3-x86_64.AppDir/opt/python3.8/lib/python3.8/site-packages/h5py/.libs/
# LIBS= $(H5PY_LIBS)/libhdf5-2d27eb21.so.103.0.0 -Wl,-rpath,$(H5PY_LIBS)

all: bin/test-hdf5 bin/test-charles

bin/test-%: test/test-%.c src/antenna.c src/antenna.h
	@mkdir -p bin
	$(CC) -o $@ $(CFLAGS) -I src $^ $(LIBS)


bin/test-charles: src/hdf5_decoder.c src/GRAND_antenna.h
	@mkdir -p bin
	$(CC) -o $@ $(CFLAGS) $< $(LIBS)

clean:
	rm -f bin
