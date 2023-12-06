#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <time.h>
#include "sine.h"
#include "cosine.h"
#include "inv_matrix.h"

// reference: [1]Labiod S, Guerra T M. Direct and indirect adaptive fuzzy control for a class of MIMO nonlinear systems[M]. INTECH Open Access Publisher, 2010.

// global variables declaration
#define PI 3.14159
#define ARRAY_SIZE 30000  // sampling times
#define n_joint 2         // number of robot manipulator joints
#define n_membership 3     // number of membership functions
#define dimension 81       // dimension of fuzzy controller parameter

static double Ts = 0.001;  // sampling period
static double t0 = 0.0;    // start time
static double t1 = 30.0;   // end time
static double l;           // length of robot manipulator joints
static double m[n_joint];  // mass of robot manipulator joints
static double kv[n_joint]; // coefficient of viscous friction
static double kc[n_joint]; // coefficient of coulomb friction
static double g = 9.41;    // gravitational acceleration
    
// calculate matrix multiplication
void matrix_Multi(double *C, double *A, double *B, int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        for (int k = 0; k < cols2; k++){
            *(C + j * cols2 + k) = 0.0;
            for (int g = 0; g < cols1; g++){
                *(C + j * cols2 + k) += *(A + j * cols1 + g) * *(B + g * cols2 + k);
            }
        }
    }
}

// calculate the transpose of matrix
void matrix_transpose(int rows, int cols, double matrix[rows][cols], double result[cols][rows]){
    for (int j = 0; j < rows; j++){
        for (int k = 0; k < cols; k++){
            result[k][j] = matrix[j][k];
        }
    }
}

// symbolic function
double sign(double a) {
    if (a > 0) {
        return 1.0;
    } else if (a < 0) {
        return -1.0;
    } else {
        return 0.0;
    }
}

struct _archive{
    double q1_archive[ARRAY_SIZE];
    double dq1_archive[ARRAY_SIZE];
    double q2_archive[ARRAY_SIZE];
    double dq2_archive[ARRAY_SIZE];
    double error1_archive[ARRAY_SIZE];
    double error2_archive[ARRAY_SIZE];
    double error1_velocity_archive[ARRAY_SIZE];
    double error2_velocity_archive[ARRAY_SIZE];
    double control1_archive[ARRAY_SIZE];
    double control2_archive[ARRAY_SIZE];
} archive;

Data q1_desired, dq1_desired, ddq1_desired, qd1_1, qd1_2, dqd1_1, dqd1_2, ddqd1_1, ddqd1_2;
Data q2_desired, dq2_desired, ddq2_desired, qd2_1, qd2_2, dqd2_1, dqd2_2, ddqd2_1, ddqd2_2;

struct Amp{
    double qd1_1, qd1_2;
    double dqd1_1, dqd1_2;
    double ddqd1_1, ddqd1_2;
    double qd2_1, qd2_2;
    double dqd2_1, dqd2_2;
    double ddqd2_1, ddqd2_2;
};

struct M0{
    double qd1_1, qd1_2;
    double dqd1_1, dqd1_2;
    double ddqd1_1, ddqd1_2;
    double qd2_1, qd2_2;
    double dqd2_1, dqd2_2;
    double ddqd2_1, ddqd2_2;
};

struct B0{
    double qd1_1, qd1_2;
    double dqd1_1, dqd1_2;
    double ddqd1_1, ddqd1_2;
    double qd2_1, qd2_2;
    double dqd2_1, dqd2_2;
    double ddqd2_1, ddqd2_2;
};

void SystemInput(Data *q1_desired, Data *qd1_1, Data *qd1_2, Data *dq1_desired, Data *dqd1_1, Data *dqd1_2, Data *ddq1_desired, Data *ddqd1_1, Data *ddqd1_2, Data *q2_desired, Data *qd2_1, Data *qd2_2, Data *dq2_desired, Data *dqd2_1, Data *dqd2_2, Data *ddq2_desired, Data *ddqd2_1, Data *ddqd2_2, double Ts, double t0, double t1){

    struct Amp amp; // amplitude
    amp.qd1_1 = PI;
    amp.qd1_2 = PI * 0.1;
    amp.dqd1_1 = PI * 0.5;
    amp.dqd1_2 = PI * 0.1 * 2;
    amp.ddqd1_1 = -PI * pow(0.5, 2);
    amp.ddqd1_2 = -PI * 0.1 * pow(2, 2);
    amp.qd2_1 = PI * 0.5;
    amp.qd2_2 = PI * 0.1;
    amp.dqd2_1 = PI * 0.5;
    amp.dqd2_2 = PI * 0.1 * 3;
    amp.ddqd2_1 = -PI * 0.5;
    amp.ddqd2_2 = -PI * 0.1 * pow(3, 2);

    struct M0 m0; // angular frequency
    m0.qd1_1 = 0.5;
    m0.qd1_2 = 2;
    m0.dqd1_1 = 0.5;
    m0.dqd1_2 = 2;
    m0.ddqd1_1 = 0.5;
    m0.ddqd1_2 = 2;
    m0.qd2_1 = 1;
    m0.qd2_2 = 3;
    m0.dqd2_1 = 1;
    m0.dqd2_2 = 3;
    m0.ddqd2_1 = 1;
    m0.ddqd2_2 = 3;

    struct B0 b0; // vertical shift
    b0.qd1_1 = 0;
    b0.qd1_2 = 0;
    b0.dqd1_1 = 0;
    b0.dqd1_2 = 0;
    b0.ddqd1_1 = 0;
    b0.ddqd1_2 = 0;
    b0.qd2_1 = 0;
    b0.qd2_2 = 0;
    b0.dqd2_1 = 0;
    b0.dqd2_2 = 0;
    b0.ddqd2_1 = 0;
    b0.ddqd2_2 = 0;

    sine(qd1_1, Ts, t0, t1, amp.qd1_1, m0.qd1_1, b0.qd1_1);
    sine(qd1_2, Ts, t0, t1, amp.qd1_2, m0.qd1_2, b0.qd1_2);
    cosine(dqd1_1, Ts, t0, t1, amp.dqd1_1, m0.dqd1_1, b0.dqd1_1);
    cosine(dqd1_2, Ts, t0, t1, amp.dqd1_2, m0.dqd1_2, b0.dqd1_2);
    sine(ddqd1_1, Ts, t0, t1, amp.ddqd1_1, m0.ddqd1_1, b0.ddqd1_1);
    sine(ddqd1_2, Ts, t0, t1, amp.ddqd1_2, m0.ddqd1_2, b0.ddqd1_2);
    sine(qd2_1, Ts, t0, t1, amp.qd2_1, m0.qd2_1, b0.qd2_1);
    sine(qd2_2, Ts, t0, t1, amp.qd2_2, m0.qd2_2, b0.qd2_2);
    cosine(dqd2_1, Ts, t0, t1, amp.dqd2_1, m0.dqd2_1, b0.dqd2_1);
    cosine(dqd2_2, Ts, t0, t1, amp.dqd2_2, m0.dqd2_2, b0.dqd2_2);
    sine(ddqd2_1, Ts, t0, t1, amp.ddqd2_1, m0.ddqd2_1, b0.ddqd2_1);
    sine(ddqd2_2, Ts, t0, t1, amp.ddqd2_2, m0.ddqd2_2, b0.ddqd2_2);

    q1_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    q2_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    dq1_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    dq2_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    ddq1_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    ddq2_desired->y = malloc(sizeof(double) * ARRAY_SIZE);

    for (int j = 0; j < ARRAY_SIZE; j++){
        q1_desired->y[j] = qd1_1->y[j] + qd1_2->y[j];       // desired angular position of joint 1
        q2_desired->y[j] = qd2_1->y[j] + qd2_2->y[j];       // desired angular position of joint 2
        dq1_desired->y[j] = dqd1_1->y[j] + dqd1_2->y[j];    // desired angular velocity of joint 1
        dq2_desired->y[j] = dqd2_1->y[j] + dqd2_2->y[j];    // desired angular velocity of joint 2
        ddq1_desired->y[j] = ddqd1_1->y[j] + ddqd1_2->y[j]; // desired angular acceleration of joint 1
        ddq2_desired->y[j] = ddqd2_1->y[j] + ddqd2_2->y[j]; // desired angular acceleration of joint 2
    }
}

struct _system_state{
    double q[n_joint];   // actual angular displacement
    double dq[n_joint];  // actual angular velocity
    double ddq[n_joint]; // actual angular acceleration
} system_state;

// linguistic variable in membership function
struct _dynamics{
    double M[n_joint][n_joint]; // inertia matrix for manipulator dynamics equation
    double C[n_joint][n_joint]; // Coriolis and centrifugal matrix for manipulator dynamics equation
    double G[n_joint];          // gravitational torque matrix for manipulator dynamics equation
    double Fv[n_joint];         // viscous friction torque
    double Fc[n_joint];         // coulomb friction torque
} dynamics;

typedef struct{
    double negative;
    double medium;
    double positive;
} membership;

// membership function of input vectors in fuzzy sets
typedef struct{
    membership error1;            // membership function of tracking error of joint 1
    membership error1_derivative; // membership function of velocity tracking error of joint 1
    membership error2;            // membership function of tracking error of joint 2
    membership error2_derivative; // membership function of velocity tracking error of joint 2
} membershipfunction;

// pointer to access membership grade
double membershipgrade(membership* member, int index) {
    switch (index) {
        case 0:
            return member->negative;
        case 1:
            return member->medium;
        case 2:
            return member->positive;
        default:
            return 0;
    }
}

// calculate Gaussian membership value
double gaussian(double x, double mean, double sigma) {
    return exp(-0.5 * pow((x - mean) / sigma, 2));
}

struct _controller{
    double controller_u[12];
    double controller_out[n_joint];
    double z[4];                                 // input variable of fuzzy controller
    double error[n_joint];                       // angular position tracking error
    double error_velocity[n_joint];              // angular velocity tracking error
    double error_second_derivative[n_joint];     // second order derivative of tracking error
    double error_state[n_joint];                 // filtered tracking error
    double error_state_derivative[n_joint];      // derivative of filtered tracking error
    double xi[dimension];                        // fuzzy basis function xi of fuzzy controller
    double xi_numerator[dimension];              // numerator of fuzzy basis function xi of fuzzy controller
    double xi_denominator;                       // denominator of fuzzy basis function xi of fuzzy controller
    double eta;                                  // learning rate
    double delta;                                // delta-modification term parameter
    double lambda[n_joint];                      // filtered tracking error parameter
    double epsilon0;                             // small positive constant
    double control[n_joint];                     // control volume
    double K[n_joint][n_joint];                  // control gain of filter tracking error
    double K0[n_joint][n_joint];                 // control gain of hyperbolic tangent function term
    double theta[dimension][n_joint];            // fuzzy controller parameter
    double theta_derivative[dimension][n_joint]; // derivative of fuzzy controller parameter
} controller;

void CONTROLLER_init(){
    system_state.q[0] = 0.1;
    system_state.dq[0] = 2.0;
    system_state.ddq[0] = 0.0;
    system_state.q[1] = 0.1;
    system_state.dq[1] = 2.5;
    system_state.ddq[1] = 0.0;
    controller.controller_u[0] = q1_desired.y[0];
    controller.controller_u[1] = dq1_desired.y[0];
    controller.controller_u[2] = ddq1_desired.y[0];
    controller.controller_u[3] = q2_desired.y[0];
    controller.controller_u[4] = dq2_desired.y[0];
    controller.controller_u[5] = ddq2_desired.y[0];
    controller.controller_u[6] = system_state.q[0];
    controller.controller_u[7] = system_state.dq[0];
    controller.controller_u[8] = system_state.ddq[0];
    controller.controller_u[9] = system_state.q[1];
    controller.controller_u[10] = system_state.dq[1];
    controller.controller_u[11] = system_state.ddq[1];

    controller.lambda[0] = 1.0; // filtered tracking error parameter
    controller.lambda[1] = 1.0;
    controller.eta = 10;        // learning rate
    controller.epsilon0 = 0.01; // small positive constant
    controller.delta = 0.001;   // delta-modification term parameter

    for (int j = 0; j < n_joint; j++) {
        for (int k = 0; k < n_joint; k++){
            if (j == k){
                controller.K[j][k] = (j == 0) ? 100.0 : 50.0;  // control gain of filter tracking error
                controller.K0[j][k] = (j == 0) ? 100.0 : 50.0; // control gain of hyperbolic tangent function term
            }
            else{
                controller.K[j][k] = 0.0;
                controller.K0[j][k] = 0.0;
            }
        }
    }

    for (int j = 0; j < n_joint; j++) {
        for (int k = 0; k < dimension; k++){
            controller.theta[k][j] = 0.0; // initialize fuzzy controller parameter
        }
    }

}

double CONTROLLER_realize(int i){
    controller.controller_u[0] = q1_desired.y[i];      // desired angular position of joint 1
    controller.controller_u[1] = dq1_desired.y[i];     // desired angular velocity of joint 1
    controller.controller_u[2] = ddq1_desired.y[i];    // desired angular acceleration of joint 1
    controller.controller_u[3] = q2_desired.y[i];      // desired angular position of joint 2
    controller.controller_u[4] = dq2_desired.y[i];     // desired angular velocity of joint 2
    controller.controller_u[5] = ddq1_desired.y[i];    // desired angular acceleration of joint 2
    controller.controller_u[6] = system_state.q[0];    // actual angular position of joint 1
    controller.controller_u[7] = system_state.dq[0];   // actual angular velocity of joint 1
    controller.controller_u[8] = system_state.ddq[0];  // actual angular acceleration of joint 1
    controller.controller_u[9] = system_state.q[1];    // actual angular position of joint 2
    controller.controller_u[10] = system_state.dq[1];  // actual angular velocity of joint 2
    controller.controller_u[11] = system_state.ddq[1]; // actual angular acceleration of joint 2

    archive.q1_archive[i] = controller.controller_u[6];
    archive.dq1_archive[i] = controller.controller_u[7];
    archive.q2_archive[i] = controller.controller_u[9];
    archive.dq2_archive[i] = controller.controller_u[10];

    controller.error[0] = q1_desired.y[i] - system_state.q[0];                       // angular position tracking error of link 1
    controller.error[1] = q2_desired.y[i] - system_state.q[1];                       // angular position tracking error of link 2
    controller.error_velocity[0] = dq1_desired.y[i] - system_state.dq[0];            // angular velocity tracking error of link 1
    controller.error_velocity[1] = dq2_desired.y[i] - system_state.dq[1];            // angular velocity tracking error of link 2
    controller.error_second_derivative[0] = ddq1_desired.y[i] - system_state.ddq[0]; // angular acceleration tracking error of link 1
    controller.error_second_derivative[1] = ddq2_desired.y[i] - system_state.ddq[1]; // angular acceleration tracking error of link 2

    archive.error1_archive[i] = controller.error[0];
    archive.error1_velocity_archive[i] = controller.error_velocity[0];
    archive.error2_archive[i] = controller.error[1];
    archive.error2_velocity_archive[i] = controller.error_velocity[1];

    // input variable of fuzzy controller
    for (int j = 0; j < 2; j++){
        controller.z[2 * j] = controller.error[j];              // tracking error of link j+1
        controller.z[2 * j + 1] = controller.error_velocity[j]; // velocity tracking error of link j+1
    }

    membershipfunction membership;

    // membership function of error1
    membership.error1.negative = gaussian(controller.z[0], -1.25, 0.6);
    membership.error1.medium = gaussian(controller.z[0], 0.0, 0.6);
    membership.error1.positive = gaussian(controller.z[0], 1.25, 0.6);

    // membership function of error1_derivative
    membership.error1_derivative.negative = gaussian(controller.z[1], -1.25, 0.6);
    membership.error1_derivative.medium = gaussian(controller.z[1], 0.0, 0.6);
    membership.error1_derivative.positive = gaussian(controller.z[1], 1.25, 0.6);

    // membership function of error2
    membership.error2.negative = gaussian(controller.z[2], -1.25, 0.6);
    membership.error2.medium = gaussian(controller.z[2], 0.0, 0.6);
    membership.error2.positive = gaussian(controller.z[2], 1.25, 0.6);

    // membership function of error2_derivative
    membership.error2_derivative.negative = gaussian(controller.z[3], -1.25, 0.6);
    membership.error2_derivative.medium = gaussian(controller.z[3], 0.0, 0.6);
    membership.error2_derivative.positive = gaussian(controller.z[3], 1.25, 0.6);

    // calculate the numerator and denominator of fuzzy basis function xi
    controller.xi_denominator = 0;
    int index = 0;
    for (int j = 0; j < n_membership; j++){
        for (int k = 0; k < n_membership; k++){
            for (int l = 0; l < n_membership; l++){
                for (int m = 0; m < n_membership; m++){
                    double mu1 = membershipgrade(&membership.error1, j);
                    double mu2 = membershipgrade(&membership.error1_derivative, k);
                    double mu3 = membershipgrade(&membership.error2, l);
                    double mu4 = membershipgrade(&membership.error2_derivative, m);
                    controller.xi_numerator[index] = mu1 * mu2 * mu3 * mu4; // numerator of fuzzy basis function xi
                    controller.xi_denominator += mu1 * mu2 * mu3 * mu4;     // denominator of fuzzy basis function xi
                    index++;
                }
            }
        }
    }

    // fuzzy basis function xi
    for (int j = 0; j < dimension; j++){
        if (controller.xi_denominator != 0) {
            controller.xi[j] = controller.xi_numerator[j] / controller.xi_denominator;
        } else {
            printf("error: denominator of xi is zero.\n");
            break;
        }
    }

    for (int j = 0; j < n_joint; j++) {
        controller.error_state[j] = controller.error[j] + controller.lambda[j] * controller.error_velocity[j]; // filtered tracking error
        controller.error_state_derivative[j] = controller.error_velocity[j] + controller.lambda[j] * controller.error_second_derivative[j]; // derivative of filtered tracking error
    }

    double K_error_state[n_joint], tanh_error_state[n_joint], K0_tanh_error_state[n_joint], temp1[n_joint], temp2[dimension][n_joint];

    // calculate control gain multiplied by filter tracking error
    matrix_Multi((double *)K_error_state, (double *)controller.K, (double *)controller.error_state, n_joint, n_joint, 1);

    // calculate hyperbolic tangent function value of ratio of filter tracking error to scaling factor epsilon
    for (int j = 0; j < n_joint; j++) {
        tanh_error_state[j] = tanh(controller.error_state[j] / controller.epsilon0);
    }

    // calculate control gain multiplied by hyperbolic tangent function value of filter tracking error
    matrix_Multi((double *)K0_tanh_error_state, (double *)controller.K0, (double *)tanh_error_state, n_joint, n_joint, 1);

    // calculate sum of filter tracking error derivative term and filter tracking error term and hyperbolic tangent function term multiplied by learning rate eta
    for (int j = 0; j < n_joint; j++) {
        temp1[j] = controller.eta * (controller.error_state_derivative[j] + K_error_state[j] + K0_tanh_error_state[j]);
    }

    // calculate eta * xi' * (error_state_derivative + K * error_state + K0 * tanh(error_state./epsilon0))'
    matrix_Multi((double *)temp2, (double *)controller.xi, (double *)temp1, dimension, 1, n_joint);

    // calculate derivative of fuzzy controller parameter theta
    for (int j = 0; j < n_joint; j++) {
        for (int k = 0; k < dimension; k++) {
            controller.theta_derivative[k][j] = temp2[k][j] - controller.eta * controller.delta * controller.theta[k][j];
        }
    }

    // update fuzzy controller parameter
    for (int j = 0; j < n_joint; j++) {
        for (int k = 0; k < dimension; k++) {
            controller.theta[k][j] += controller.theta_derivative[k][j] * Ts;
        }
    }

    // control amount is equal to fuzzy basis function xi multiplied by fuzzy controller parameter theta
    matrix_Multi((double *)controller.control, (double *)controller.xi, (double *)controller.theta, 1, dimension, n_joint); // fuzzy controller output

    controller.controller_out[0] = controller.control[0];
    controller.controller_out[1] = controller.control[1];
    archive.control1_archive[i] = controller.control[0];
    archive.control2_archive[i] = controller.control[1];
}

struct _plant{
    double plant_u[2];
    double plant_out[6];
} plant;

void PLANT_init(){
    system_state.q[0] = 0.1;
    system_state.dq[0] = 2.0;
    system_state.ddq[0] = 0.0;
    system_state.q[1] = 0.1;
    system_state.dq[1] = 2.5;
    system_state.ddq[1] = 0.0;
    plant.plant_out[0] = system_state.q[0];
    plant.plant_out[1] = system_state.dq[0];
    plant.plant_out[2] = system_state.ddq[0];
    plant.plant_out[3] = system_state.q[1];
    plant.plant_out[4] = system_state.dq[1];
    plant.plant_out[5] = system_state.ddq[1];
    l = 1.0;     // length of linl
    m[0] = 1.0;  // mass of link 1
    m[1] = 2.0;  // mass of link 2
    kv[0] = 0.3; // coefficient of viscous friction
    kv[1] = 0.5;
    kc[0] = 0.2; // coefficient of coulomb friction
    kc[1] = 0.5;
}

double PLANT_realize(int i){
    plant.plant_u[0] = controller.controller_out[0];
    plant.plant_u[1] = controller.controller_out[1];
    // inertia matrix for manipulator dynamics equation
    dynamics.M[0][0] = 1.0/3.0 * m[0] * pow(l, 2) + 4.0/3.0 * m[1] * pow(l, 2) + m[1] * pow(l, 2) * cos(system_state.q[1]);
    dynamics.M[0][1] = 1.0/3.0 * m[1] * pow(l, 2) + 1.0/2.0 * m[1] * pow(l, 2) * cos(system_state.q[1]);
    dynamics.M[1][0] = dynamics.M[0][1];
    dynamics.M[1][1] = 1.0/3.0 * m[1] * pow(l, 2);

    // Coriolis and centrifugal matrix for manipulator dynamics equation 
    dynamics.C[0][0] = -m[1] * pow(l, 2) * sin(system_state.q[1]) * system_state.dq[1];
    dynamics.C[0][1] = -1.0/2.0 * m[1] * pow(l, 2) * sin(system_state.q[1]) * system_state.dq[1];
    dynamics.C[1][0] = 1.0/2.0 * m[1] * pow(l, 2) * sin(system_state.q[1]) * system_state.dq[0];
    dynamics.C[1][1] = 0.0;

    // gravitational torque matrix for manipulator dynamics equation
    dynamics.G[0] = 1.0/2.0 * m[0] * g * l * cos(system_state.q[0]) + 1.0/2.0 * m[1] * g * l * cos(system_state.q[0] + system_state.q[1]) + m[1] * g * l * cos(system_state.q[0]);
    dynamics.G[1] = 1.0/2.0 * m[1] * g * l * cos(system_state.q[0] + system_state.q[1]);

    // viscous friction torque
    dynamics.Fv[0] = kv[0] * system_state.dq[0];
    dynamics.Fv[1] = kv[1] * system_state.dq[1];

    // coulomb friction torque
    dynamics.Fc[0] = kc[0] * sign(system_state.dq[0]);
    dynamics.Fc[1] = kc[1] * sign(system_state.dq[1]);

    // for (int j = 0; j < 2; j++) {
    //     printf("dynamics.Fc[%d]: %lf\n", j, dynamics.Fc[j]);
    // }

    double inv_M[n_joint][n_joint], C_dq[n_joint], control_Cdq_G_Fv_Fc[n_joint];
    inv_matrix(inv_M, dynamics.M, n_joint); // calculate inverse of inertia matrix for manipulator dynamics equation
    matrix_Multi((double *)C_dq, (double *)dynamics.C, (double *)system_state.dq, 2, 2, 1); // calculate C multiplied by actual angular velocity

    for (int j = 0; j < n_joint; j++){
        control_Cdq_G_Fv_Fc[j] = controller.control[j] - C_dq[j] - dynamics.G[j] - dynamics.Fv[j] - dynamics.Fc[j];
    }

    // actual manipulator dynamics system
    matrix_Multi((double *)system_state.ddq, (double *)inv_M, (double *)control_Cdq_G_Fv_Fc, 2, 2, 1);

    system_state.dq[0] = system_state.dq[0] + system_state.ddq[0] * Ts;
    system_state.dq[1] = system_state.dq[1] + system_state.ddq[1] * Ts;
    system_state.q[0] = system_state.q[0] + system_state.dq[0] * Ts;
    system_state.q[1] = system_state.q[1] + system_state.dq[1] * Ts;

    plant.plant_out[0] = system_state.q[0];
    plant.plant_out[1] = system_state.dq[0];
    plant.plant_out[2] = system_state.ddq[0];
    plant.plant_out[3] = system_state.q[1];
    plant.plant_out[4] = system_state.dq[1];
    plant.plant_out[5] = system_state.ddq[1];

}

void saveArchiveToTxt(double *archive, int size, const char *filename){

    FILE *file = fopen(filename, "w+");

    if (file == NULL){
        perror("Failed to open file");
        exit(1);
    }
    else{
        for (int i = 0; i < size; i++){
            fprintf(file, "%lf\n", archive[i]);
        }
        fclose(file);
        printf("Saved to file %s\n", filename);
    }
}

void saveArchive(){

    saveArchiveToTxt(q1_desired.y, ARRAY_SIZE, "../report/qd1.txt");
    saveArchiveToTxt(dq1_desired.y, ARRAY_SIZE, "../report/dqd1.txt");
    saveArchiveToTxt(archive.q1_archive, ARRAY_SIZE, "../report/q1.txt");
    saveArchiveToTxt(archive.dq1_archive, ARRAY_SIZE, "../report/dq1.txt");
    saveArchiveToTxt(q2_desired.y, ARRAY_SIZE, "../report/qd2.txt");
    saveArchiveToTxt(dq2_desired.y, ARRAY_SIZE, "../report/dqd2.txt");
    saveArchiveToTxt(archive.q2_archive, ARRAY_SIZE, "../report/q2.txt");
    saveArchiveToTxt(archive.dq2_archive, ARRAY_SIZE, "../report/dq2.txt");
    saveArchiveToTxt(archive.error1_archive, ARRAY_SIZE, "../report/error1.txt");
    saveArchiveToTxt(archive.error1_velocity_archive, ARRAY_SIZE, "../report/error1_velocity.txt");
    saveArchiveToTxt(archive.error2_archive, ARRAY_SIZE, "../report/error2.txt");
    saveArchiveToTxt(archive.error2_velocity_archive, ARRAY_SIZE, "../report/error2_velocity.txt");
    saveArchiveToTxt(archive.control1_archive, ARRAY_SIZE, "../report/control1.txt");
    saveArchiveToTxt(archive.control2_archive, ARRAY_SIZE, "../report/control2.txt");
}

int main(){


    SystemInput(&q1_desired, &qd1_1, &qd1_2, &dq1_desired, &dqd1_1, &dqd1_2, &ddq1_desired, &ddqd1_1, &ddqd1_2, &q2_desired, &qd2_1, &qd2_2, &dq2_desired, &dqd2_1, &dqd2_2, &ddq2_desired, &ddqd2_1, &ddqd2_2, Ts, t0, t1);
    CONTROLLER_init(); // initialize controller parameter
    PLANT_init();      // initialize plant parameter

    for (int i = 0; i < ARRAY_SIZE; i++){
    // for (int i = 0; i < 5; i++){
        double time = i * Ts + t0;
        printf("time at step %d: %f\n", i, time);
        CONTROLLER_realize(i);
        PLANT_realize(i);
    }

    saveArchive();

    return 0;
}
