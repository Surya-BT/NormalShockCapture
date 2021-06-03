# NormalShockCapture
In this project, I used simulated the normal shock capture in a converging-diverging nozzle (quasi 1D) using the MacCorMack Method. 

The artificial viscosity or numerical dissipiation affect the accuracy of the solution, but also help in stabilizing it. These become particularly useful when capturing flows with very large gradients such as shock waves.

I used a 4th order numerical dissipiation term to capture the shock waves and compared with the solution without numerical dissipiation.
It can be observed from the video that the addition of artificial viscosity helped remove the spurious waves near the shock and stabilize the solution.

It should also be noted that the usage of large amounts of numerical dissipiation can smooth the gradients and cost the accuracy of the solution. Hence the selection of the proportion of artificial viscosity is always dependent on the problem-type.

The analytical solutions were obtained from the book referenced below

References:
Computational Fluid Dynamics: The basics with applications by JD Anderson
Link: Anderson, J.. “Computational fluid dynamics : the basics with applications.” (1995).
