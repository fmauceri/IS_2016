# IS_2016
###### Grace Handler, Stephen Majercik, Frank Mauceri, Jack Truskowski


![Boids](https://cloud.githubusercontent.com/assets/11000833/12837197/170d39a2-cb8e-11e5-84c1-273e1cb236c9.png)


Instructions for running the simulation can be found in the `ControlPanel...maxpat` file

####Future Work:
- Pattr: how to save, finding cool combos

- Visually connect boids in same neighborhood (output a matrix from external with each boid and its neighbors) This will make neighborhoods / clusters more visually obvious
- Revisit the statistics used to characterize flock behavior  
- Stipulate where boids are born in space when they are added 

###Questions:
- Amount of attraction based on attractor, not flocks? Or both?
- For clustering: what if we passed a matrix with planes per boid and within plane first line is xyz of that boid and following lines xyz of boids in neighborhood (PROBLEM: you'd be looking at every twice)
- Revisit the statistics used to characterize flock behavior
- Stipulate where boids are born in space when they are added

- NOTE: repel and edge parameters are currently not hooked up
- Theres a weird bug that only occurs sometimes, if boids are in all flocks, where they're not drawn (I think)



