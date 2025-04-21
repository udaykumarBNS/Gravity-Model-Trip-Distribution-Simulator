#include <bits/stdc++.h>
using namespace std;

const double range = 0.001; // Convergence threshold
const double ALPHA = 2.0;   // Given alpha value
const double EPSILON = 1e-6; // Small value to avoid division by zero
const int MAX_ITER = 1000;  // Maximum iterations for convergence

// Function for Production Constrained Gravity Model
void productionConstrainedModel(int n, vector<double> &P, vector<vector<double>> &d) {
    vector<double> Ai(n, 1.0);
    vector<vector<double>> T(n, vector<double>(n, 0.0));

    // Compute Ai values
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += P[j] * pow(d[i][j], ALPHA);
        }
        Ai[i] = (sum > EPSILON) ? 1.0 / sum : 1.0;
    }

    // Compute trip distribution matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (d[i][j] < EPSILON) {
                T[i][j] = 0.0;
            } else {
                T[i][j] = (Ai[i] * P[i] * P[j]) / pow(d[i][j], ALPHA);
            }
        }
    }

    // Display the computed OD matrix
    cout << "\nOD Trip Distribution Matrix (T_ij) - Production Constrained:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << fixed << setprecision(2) << T[i][j] << " ";
        }
        cout << endl;
    }
}

// Function for Attraction Constrained Gravity Model
void attractionConstrainedModel(int n, vector<double> &A, vector<vector<double>> &d) {
    vector<double> Bj(n, 1.0);
    vector<vector<double>> T(n, vector<double>(n, 0.0));

    // Compute Bj values
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += A[i] * pow(d[i][j], ALPHA);
        }
        Bj[j] = (sum > EPSILON) ? 1.0 / sum : 1.0;
    }

    // Compute trip distribution matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (d[i][j] < EPSILON) {
                T[i][j] = 0.0;
            } else {
                T[i][j] = (Bj[j] * A[i] * A[j]) / pow(d[i][j], ALPHA);
            }
        }
    }

    // Display the computed OD matrix
    cout << "\nOD Trip Distribution Matrix (T_ij) - Attraction Constrained:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << fixed << setprecision(2) << T[i][j] << " ";
        }
        cout << endl;
    }
}

// Function for Doubly Constrained Gravity Model
void gravityModel(int n, vector<double> &P, vector<double> &A, vector<vector<double>> &d) {
    vector<double> Ai(n, 1.0);
    vector<double> Bj(n, 1.0);
    vector<vector<double>> T(n, vector<double>(n, 0.0));

    bool converged = false;
    vector<double> Ai_prev(n, 0.0);
    vector<double> Bj_prev(n, 0.0);
    int iter_count = 0;

    // Ensure P and A are nonzero to prevent NaN issues
    for (int i = 0; i < n; i++) {
        if (P[i] < EPSILON) P[i] = EPSILON;
        if (A[i] < EPSILON) A[i] = EPSILON;
    }

    // Initial better guess for Bj
    for (int j = 0; j < n; j++) {
        Bj[j] = 1.0 / (A[j] + EPSILON);
    }

    while (!converged && iter_count < MAX_ITER) {
        iter_count++;

        // Compute Ai values
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += Bj[j] * A[j] / pow(max(d[i][j], EPSILON), ALPHA);
            }
            Ai[i] = (sum > EPSILON) ? 1.0 / sum : 1.0;
        }

        // Compute Bj values
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += Ai[i] * P[i] / pow(max(d[i][j], EPSILON), ALPHA);
            }
            Bj[j] = (sum > EPSILON) ? 1.0 / sum : 1.0;
        }

        // Check for convergence
        converged = true;
        for (int i = 0; i < n; i++) {
            if (fabs(Ai[i] - Ai_prev[i]) > range || fabs(Bj[i] - Bj_prev[i]) > range) {
                converged = false;
                break;
            }
        }

        Ai_prev = Ai;
        Bj_prev = Bj;
    }

    // Compute trip distribution matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (d[i][j] < EPSILON) {
                T[i][j] = 0.0;
            } else {
                T[i][j] = (Ai[i] * Bj[j] * P[i] * A[j]) / pow(max(d[i][j], EPSILON), ALPHA);
            }
        }
    }

    // Display the computed OD matrix
    cout << "\nOD Trip Distribution Matrix (T_ij) - Doubly Constrained:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << fixed << setprecision(2) << ceil(T[i][j]) << " ";
        }
        cout << endl;
    }
}

int main() {
    int choice, n;
    cout << "Choose the model:\n";
    cout << "1 - Production Constrained Gravity Model\n";
    cout << "2 - Attraction Constrained Gravity Model\n";
    cout << "3 - Doubly Constrained Gravity Model\n";
    cout<<"Enter the model: "<<"";
    cin >> choice;

    cout << "Enter the number of zones: ";
    cin >> n;

    vector<double> P(n), A(n);
    vector<vector<double>> d(n, vector<double>(n));

    if (choice == 1) {
        cout << "Enter production values for each zone:\n";
        for (int i = 0; i < n; i++) {
            cout << "Enter P[" << i + 1 << "]: ";
            cin >> P[i];
        }

        cout << "Enter the distance matrix:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << "Enter d[" << i + 1 << "][" << j + 1 << "]: ";
                cin >> d[i][j];
            }
        }

        productionConstrainedModel(n, P, d);

    } else if (choice == 2) {
        cout << "Enter attraction values for each zone:\n";
        for (int j = 0; j < n; j++) {
            cout << "Enter A[" << j + 1 << "]: ";
            cin >> A[j];
        }

        cout << "Enter the distance matrix:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << "Enter d[" << i + 1 << "][" << j + 1 << "]: ";
                cin >> d[i][j];
            }
        }

        attractionConstrainedModel(n, A, d);

    } else if (choice == 3) {
        cout << "Enter production values for each zone:\n";
        for (int i = 0; i < n; i++) {
            cout << "Enter P[" << i + 1 << "]: ";
            cin >> P[i];
        }

        cout << "Enter attraction values for each zone:\n";
        for (int i = 0; i < n; i++) {
            cout << "Enter A[" << i + 1 << "]: ";
            cin >> A[i];
        }

        cout << "Enter the distance matrix:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << "Enter d[" << i + 1 << "][" << j + 1 << "]: ";
                cin >> d[i][j];
            }
        }

        gravityModel(n, P, A, d);

    } else {
        cout << "Invalid choice!" << endl;
    }

    return 0;
}