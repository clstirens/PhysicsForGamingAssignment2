// PhysicsProject2.cpp : Defines the entry point for the console application. 
// Developed by Charles Stirens for Physics Project 2, 2d object flying through resistive medium
// Program developed to accept changes to firing angle, initial velocity, time, and terminal velocity found at the top of the main method.
// Version 1.0 Nov 30 2016

#include "stdafx.h" 
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

#define PI 3.14159265 //Define pi to convert from degrees to radians. Note the accuracy of the value PI, this is addressed in our degreeToRad method.

using namespace std;
double degreeToRad(double param); //init degree conversion method

int main()
{
	double theta; //initialize theta, angle in degrees in which the projectile is fired.
	string dim1 = "output1d.csv";
	string dim2 = "output2d.csv";
	string dim2drag = "output2dDrag.csv";
	double firingAngle = 50; // create value for angle theta. Change this to alter the test case.
	double v0 = 60; //meters per second that our projectile is initially launched, change this to alter test case.
	double x0 = 0; // initial value of x in meters. This value is used because the object has not moved yet.
	double y0 = 0; // initial value of y in meters. This value is used because the object has just been launched
	double vTerm = 60; // our designated terminal velocity in meters per second.
	double g = 9.8; // gravity = 9.8 m/s
	double t = 26; // Time of terminal velocity for our object
	double stepsize = .25; // Stepsize for accuracy, Iterate t by this. Many equations will scale based on this value.
	double dr; //initiliaze drag value to be used later.

	ofstream myfile; // inititalize to file operations

	for (int j = 0; j < 3; j++) {

		if (j == 0) //if on first test case, angle 90 is used to simulate 1d
			theta = 90; //set theta to 90 to simulate 1d object
		else
			theta = firingAngle; //sets theta to angle of choice in order to test 2d object with and without drag.
		
		double v = v0; // initialize and set v to v0 value. This will be used in the first step in the euler method

		double vx = v0*cos(degreeToRad(theta)); //find initial velocity in x by using cos of our angle times initial velocity
		double vy = v0*sin(degreeToRad(theta)); //find initial velocity in y by using sin of our angle times initial velocity

		double x = x0; // initializes and sets displacement value in x based on initial value
		double y = y0; // initalizes and sets displacement value in y based on initial value

		if (j == 1) // if on second test case, drag is set to zero
			dr = 0;
		else
			dr = (g / (vTerm * vTerm)); // b/m value, dragCoefficient, quadratic drag.

		double ax = - (dr * v * vx); //find initial value acceleration in x
		double ay = - g - (dr * v * vy); //find initial value acceleration in y

		
		/*This section changes the file we open for each test case. Each file name is tied to a string.
		These are called dim1, dim2, dim3. The top of the program contains these file names, which can be changed if necessary*/
			if (j == 0)
				myfile.open(dim1); 
			else if (j == 1)
				myfile.open(dim2);
			else
				myfile.open(dim2drag);

		myfile << "Time (t),Velocity X (Vx),Velocity Y (Vy),Displacement X,Displacement Y\n"; // table labels for excel doc
		cout << dr << "\n";

		for (int i = 0; i <= (t / stepsize); i++) { //for number of iterations based on step size. t/step in order to scale our simulation dynamically, in case another user wants to test other stepsizes.

			cout << i*stepsize << "   " << vx << "   " << vy << "   " << x << "   " << y << "\n"; //cout for debug, comment out if needed(left uncommented to ensure program finishes)
			myfile << i*stepsize << "," << vx <<  "," << vy << "," << x << ","  << y <<  "\n"; //output each variable to comma seperated file. base case uses starting values
			if (y <= 0 && i>1) //if y falls below zero AND we're not on the first step, stop processing. In later versions, make this toggle-able.
				i = ((t / stepsize) + 1);
			vx = vx + (ax * (stepsize)); //Determine (n+1) x velocity value based on previous value + acceleration in x * step size
			vy = vy + (ay * (stepsize)); //Determine (n+1) y velocity value based on previous value + acceleration in y * step size

			x = x + (vx * (stepsize)); //displacement in x = previous x displacement + (x velocity * step)
			y = y + (vy * (stepsize)); //displacement in y = previous y displacement + (y velocity * step)

			v = sqrt((vx*vx) + (vy*vy)); // find velocity based on the sqrt of x velocity^2 + y velocity^2
		
			ax = - (dr * v * vx); //acceleration in x affected by drag * velocity * x velocity
			ay = -g - (dr * v * vy); //acceleration in x affected by drag * velocity * y velocity
			
		}
	
	

		myfile.close(); //closes and saves file for use later
	}
	system("pause"); //pauses runtime window and prompts user
	return 0; //end
}

/*
This is a method created to evaluate our angle in degrees and convert to radians. 
Note that our declaration of PI at the top is only a small part of PI, and can lead to a small accuracy loss but can benefit performance.
*/
double degreeToRad(double param) { //Method for converting degrees to radians, takes parameter from program
	double result; // initialize result variable
	result = (param * PI) / 180.0; // Parameter * PI / 180 degrees. 
	return result; //Return radians value
}