# Distributed Systems for Measurements and Automation 19/20

## Project by Rubio Gomez Alvaro & Boffo Marco
Hybrid algorithm for distributing groups of robots. 
Our approach allows to distribute them homogeneously on an area characterised by a variable density function. 
The focus is on autonomous agents either able to perform sensing tasks or able to receive a desired configuration.

## Simulation

 - Open ```Simulation.m``` with MatLab

 - Set the initial parameters as desired:
    ```
    % Number of robots
    N = 9; 
    % Dimension of the picture square to cover
    dim = 100; 
    % Polygon area definition (square)
    crs = [0,0;0,dim;dim,dim;dim,0]; 
    % Radius of the robots
    R = 25 + rand(N,1)*10; 

    % Number of maximum executions of the algorithm
    iterations = 10000; 
    % May include the name of the file inside images/ folder (equal to end_image to stop after a single image)
    image_i = 0; % Incremental to change image
    end_image = 1; % Set number for final image +1 (>1 to test "moving images")
    step = 0.7; % Maximum distance that a robot can move in a direction, 0.7 so to have maximum length of 1

    % Initial position of the robots
    abs_value_x = dim * rand(N,1)
    abs_value_y = dim * rand(N,1)
    ```
