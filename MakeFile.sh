DELPHES_DIR=/home/john/Programs/Delphes-3.3.2

EXROOT_DIR=/external/ExRootAnalysis
FJ_DIR=/external

SO_FILE=Delphes

DELPHES_LIB="-Wl,-rpath,$DELPHES_DIR -L$DELPHES_DIR -lDelphes"

CUSTOM_ANALYSIS="examples/Analysis.cpp examples/processes.cpp"

CXXFLAGS="-I/home/john/Programs/root/include -I$DELPHES_DIR -I$DELPHES_DIR/classes -I$DELPHES_DIR/$EXROOT_DIR -I$DELPHES_DIR/$FJ_DIR"
LDFLAGS="-L/home/john/Programs/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lEG $DELPHES_LIB"

g++ $CXXFLAGS examples/jetMass_s.cpp $LDFLAGS $CUSTOM_ANALYSIS -o jetMass_s.exe
