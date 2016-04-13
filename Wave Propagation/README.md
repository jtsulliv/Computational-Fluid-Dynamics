# Wave Propagation

In this project three different finite difference schemes are used to discretize the governing equation for the convection-diffusion partial differential equation (PDE).  The schemes implemented include a first-order forward-time and second-order central space (FTCS) scheme for convection and diffusion, first order upwind for convection and FTCS for diffusion, and a MacCormack method for convection and second order central space for diffusion.  The convection portion of the governing equation is of the form of a hyperbolic PDE, while the diffusion term is of the form of a parabolic PDE.  

To verify the accuracy of the numerical scheme developed, a mesh refinement study is conducted.  Both the truncation error and stability conditions for the FTCS scheme are investigated.  In order to validate the solution, the results obtained from the numerical scheme are compared to those obtained from the analytical solution. 
