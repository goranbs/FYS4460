#include <unittest++/UnitTest++.h>
#include <atom.h>
#include <unitconverter.h>

TEST(Atom) {
    Atom argon({1,2,3},{1,2,3},{0.1,0.1,0.1},0.0);
    vector <double> i_pos = argon.return_initial_position();
    double truevalue = false;
    if (i_pos[0] == 1){
        if (i_pos[1] == 2){
            if (i_pos[2] == 3){
                truevalue = true;
            }
        }
    }
    CHECK(truevalue == true);
}


int main() {

    return UnitTest::RunAllTests();
}
