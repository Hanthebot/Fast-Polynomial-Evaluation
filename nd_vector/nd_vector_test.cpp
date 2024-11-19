#include "nd_vector.h"
#include <vector>

using namespace std;

enum {QUIT, EDIT, PRINT} options;

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
    
    int command, ind1, ind2, ind3;
    do {
        cout << "Choose action: " << endl
            << "0. exit" << endl
            << "1. edit value" << endl
            << "2. print" << endl;
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
                if (ind3 == -1) {
                    cout << vecs[ind1][ind2] << endl;
                    break;
                }
                cout << vecs[ind1][ind2][ind3] << endl;
                break;
        }
    } while (command);
    return 0;
}