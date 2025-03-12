#include <stdio.h>
#include <math.h>
//Define the Rosenbrock function as stated:
double rosenbrock(double x0, double x1) {
    return pow(1 - x0, 2) + 100 * pow(x1 - pow(x0, 2), 2);
}
//Define the function used to write to files
void write_to_file(double x1) { 
    FILE *datafile;
    datafile = fopen("values.txt", "w");
    //Stating the limits. The divisions of 0.04 ensuring 100 values outputted into file.
    for (double x0 = -2; x0 <= 2; x0 += 0.038) { 
        double result = rosenbrock(x0, x1);
        fprintf(datafile, "%.4lf, %lf\n", x0, result);//Format of calcualted values.
    }
    fclose(datafile);
    printf("\nValues written to file:\n\n");
}

//Function for finding the order of indicies based on their values
//Comparison checks are done to find the order of the vertex points from highest to lowest of all 3 points 
void ordering(double x1[3], int *highest, int *second_highest, int *lowest) {
    *highest = 0;
    *second_highest = 0;
    *lowest = 0;
    for (int i = 0; i < 3; i++) {
        if (x1[i] > x1[*highest]) {
            *second_highest = *highest;
            *highest = i;
        } 
        else if (x1[i] > x1[*second_highest]) {
            *second_highest = i;
        }
        if (x1[i] < x1[*lowest]) {
            *lowest = i;
        }
    }
}
//Calculating the centroid from the initial points.
//Follows the method outlined 
void calculate_centroid(double points[3][2], int highest_index, double *centroid_x,
    double *centroid_y) {
    double x_sum = 0.0, y_sum = 0.0;
    for (int i = 0; i < 3; i++) {
        if (i != highest_index) {
            x_sum += points[i][0];
            y_sum += points[i][1];
        }
    }
    *centroid_x = x_sum / 2;
    *centroid_y = y_sum / 2;
}
//Function used for reflecting when criteria for reflection met
void reflection(double points[3][2], int highest_index, double centroid_x,
    double centroid_y,
             double *reflected_x, double *reflected_y, double *func_evaluate) {
    *reflected_x = centroid_x + (centroid_x - points[highest_index][0]);
    *reflected_y = centroid_y + (centroid_y - points[highest_index][1]);
    *func_evaluate = rosenbrock(*reflected_x, *reflected_y);
}
//Function used for expanding when criteria met
void expand(double points[3][2], int highest_index, double centroid_x, double centroid_y,
            double reflected_x, double reflected_y, double func_evaluate, double *expanded_x,
            double *expanded_y, double *func_alternate, double *x1) {
    *expanded_x = centroid_x + 2.0 * (reflected_x - centroid_x);
    *expanded_y = centroid_y + 2.0 * (reflected_y - centroid_y);
    *func_alternate = rosenbrock(*expanded_x, *expanded_y);
    if (*func_alternate < func_evaluate) {
        points[highest_index][0] = *expanded_x;
        points[highest_index][1] = *expanded_y;
        *x1 = *func_alternate;
    } else {
        points[highest_index][0] = reflected_x;
        points[highest_index][1] = reflected_y;
        *x1 = func_evaluate;
    }
}
//Function used for contracting when criteria met 
void contract(double points[3][2], int highest_index, double centroid_x, double centroid_y,
              double *contracted_x, double *contracted_y, double *func_alternate, double *x1) {
    *contracted_x = centroid_x + 0.5 * (points[highest_index][0] - centroid_x);
    *contracted_y = centroid_y + 0.5 * (points[highest_index][1] - centroid_y);
    *func_alternate = rosenbrock(*contracted_x, *contracted_y);
    if (*func_alternate < *x1) {
        points[highest_index][0] = *contracted_x;
        points[highest_index][1] = *contracted_y;
        *x1 = *func_alternate;
    }
  }
//The function where criteria is tested and compared
//All the operations are also undertaken in this
void simplex_method(double points[3][2]){
    double x1[3];
    //These values are used in the movements instead of using fractions, 0.5 is used in multiplacation
    int i, j, l, highest_index, s, m, num_iterations = 0;
    double sum, centroid_x, centroid_y, reflected_x, reflected_y, func_evaluate,
        contracted_x, expanded_y, func_alternate,expanded_x,contracted_y;
    double tolerance = 1e-8; // set the tolerance value here
    // Evaluate the Rosenbrock function at the initial points
    for (i = 0; i < 3; i++) {
        x1[i] = rosenbrock(points[i][0], points[i][1]);
    }
    //Ensuring it runs until max iterations is reached
    while (num_iterations < 1000) {
        // Find the order of indices.
        ordering(x1, &highest_index, &s, &l);
        //Check to see if tolerance criteria has been reached
        if (sqrt((pow((x1[highest_index] - x1[l]),2)/2)) < tolerance) {
            break;
        }
        //Otherwise continue running through iterations
        //Calculate the centroid of the non-highest point 
        calculate_centroid(points, highest_index, &centroid_x, &centroid_y);
        // Reflect the highest point through the centroid
        reflection(points, highest_index, centroid_x, centroid_y, &reflected_x, &reflected_y, &func_evaluate);
        if (func_evaluate < x1[l]) {
            // Expand the point
            expand(points, highest_index, centroid_x, centroid_y, reflected_x, reflected_y, func_evaluate, &expanded_x, &expanded_y, &func_alternate, &x1[highest_index]);
        } else if (func_evaluate >= x1[s]) {
        // If the reflected point is worse than the second highest point, contract the simplex
            if (func_evaluate < x1[highest_index]) {
            // If the reflected point is still better than the highest point, accept it
                points[highest_index][0] = reflected_x;
                points[highest_index][1] = reflected_y;
                x1[highest_index] = func_evaluate;
            }
            // Calculate the point to which the highest point is to be contracted
            contract(points, highest_index, centroid_x, centroid_y, &contracted_x, &contracted_y, &func_alternate, &x1[highest_index]);

        } else {
        // If the reflected point is neither better nor worse than the second highest point, accept it
        points[highest_index][0] = reflected_x;
        points[highest_index][1] = reflected_y;
        x1[highest_index] = func_evaluate;
    }
    //Increment the iteration counter by 1
    num_iterations++;
}

// Output the information of interest
printf("Final Coordinates:\n(%f,%f)\n(%f,%f)\n(%f,%f)\n",points[0][0],points[0][1],points[1][0],points[1][1],points[0][2],points[2][0]);
printf("Minimum Point found at (x0,x1) = (%f, %f)\n", points[0][0], points[0][1]);
printf("Function Value, F(x0,x1), at minimum = %f\n", x1[l]);
printf("Number of iterations completed: %d\n", num_iterations);
}   

//Main function for running code and implementing functions
int main(){
    //Calling function and stating x1 is equal to 1
    write_to_file(1);
    //Defining the starting vertices of the simplex
    double points[3][2] = {{0.0, 0.0}, {2.0, 0.0}, {0.0, 2.0}};
    //Run the function for the starting coordinates specified
    simplex_method(points);
    return 0;
}
