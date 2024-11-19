#include "nd_vector.h"
#include <vector>

using namespace std;

enum {QUIT, EDIT, PRINT, SET_RANGE} options;

int main() {
    vector<int> vec;
    vec.resize(60, 0);
    span<int> span_vec{vec};

    for (size_t i = 0; i < vec.size(); ++i)
        vec[i] = i;

    size_t dim = 3;
    size_t shape_arr[] = {4, 3, 5};
    span<size_t> shape {shape_arr};
    nd_vector<int> vecs(dim, shape, span_vec);
    
    nd_vector<int> vecs2 = vecs[2][1];

    nd_vector<int> vecs3 = vecs[1];

    nd_vector<int> vvecs[] = {vecs, vecs2, vecs3};

    int command, ind1, ind2, ind3;
    do {
        cout << "Choose action: " << endl
            << "0. exit" << endl
            << "1. edit value" << endl
            << "2. print" << endl
            << "3. set range" << endl;
        cin >> command;
        switch (command) {
            case EDIT:
                cout << "Enter indices" << endl;
                cin >> ind1 >> ind2 >> ind3;
                cout << "Enter value" << endl;
                cin >> vecs[ind1][ind2][ind3].get(0);
                break;
            case PRINT:
                cout << "Enter indices to print, -1 to print" << endl;
                cin >> ind1;
                if (ind1 == -1) {
                    cout << vecs << endl;
                    break;
                }
                cin >> ind2;
                if (ind2 == -1) {
                    cout << vecs[ind1] << endl;
                    break;
                }
                cin >> ind3;
                cout << vecs[ind1][ind2][ind3] << endl;
                break;
            case SET_RANGE:
                cout << "Enter indices to set range" << endl;
                cin >> ind1;
                cout << "Enter indices of vector to copy from" << endl;
                cin >> ind2;
                vvecs[ind1].set_range(vvecs[ind2]);
                break;
        }
    } while (command);
    return 0;
}