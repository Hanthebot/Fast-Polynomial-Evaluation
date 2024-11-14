#include <iostream>
#include <span>
#include <vector>

using namespace std;

enum {QUIT, EDIT, PRINT} option;

int main() {
    vector<int> vec;
    vec.resize(15, 0);
    span<int> span_vec{vec};

    for (int i = 0; i < vec.size(); ++i)
        vec[i] = i;

    span<int> vec1 {vec.begin(), vec.begin() + 10 + 1};
    span<int> vec2 {vec.begin()+5, vec.begin() + 15};
    span<int> vecs[] = {span_vec, vec1, vec2};
    int command, option, ind;
    do {
        cout << "Choose action: " << endl
            << "0. exit" << endl
            << "1. edit value" << endl
            << "2. print" << endl;
        cin >> command;
        switch (command) {
            case EDIT:
                cout << "Enter vector (0, 1, 2) and index" << endl;
                cin >> option >> ind;
                cout << "Enter value" << endl;
                cin >> vecs[option][ind];
                break;
            case PRINT:
                cout << "Enter vector (0, 1, 2) and index" << endl;
                cin >> option;
                cout << "[ ";
                for (int v_i : vecs[option])
                    cout << v_i << " ";
                cout << "]" << endl << endl;
                break;
        }
    } while (command);
    return 0;
}