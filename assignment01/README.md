---
title: "Assignment 01"
author: "Dennys Huber"
date: \today
subject: "PAMPI-WS25/26"
header-left: "PAMPI-WS25/26"
header-right: \theauthor
footer-left: \today
footer-center: 
titlepage: false
toc: false
header-includes: |
  \usepackage{caption}
  \captionsetup[figure]{
    justification=centering,
    font=small,
    labelfont=bf
  }
---

# Assignment 1

## Part 1
Solution provided in code and can be found in the zip file.

## Part 2
The appropriate points are circled in red, representing the inter-node scaling at 20, 40, 60, and 80 processes. The other points from 1-20 processes are not appropriate for Amdahl's law fitting because within a socket all cores are competing for shared memory bandwidth, leading to saturation effects that create a bottleneck. This violates the model's assumption of scalable resources (Lecture 1, Slide 40). In order to fit Amdahl's law, we use inter-node scaling where intra-socket communication overhead is not the limiting factor. We define the first node (20 processes) as the scaling baseline, following the methodology shown in Lecture 1, Slide 40

![Appropriate data points shown in red circle](images/ambdahl.png){ width=80% }

### Serial fraction
Using the first node as our scaling baseline, we divide all speedup values by the speedup at 20 processes, resulting in the following rebased values. (Note: We now express parallelism in terms of nodes $N$, where 1 node equals 20 processes)

\begin{align}
S(1) &= \frac{12}{12} = 1\\
S(2) &= \frac{24}{12} = 2\\
S(3) &= \frac{34}{12} = \frac{17}{6}\\
S(4) &= \frac{46}{12} = \frac{23}{6}\\
\end{align}

We can now compute the serial fraction for $N = 4$ Nodes:

\begin{align}
S(N) &= \frac{1}{s + \frac{1-s}{N}} \\
\frac{23}{6} &= \frac{1}{s + \frac{1-s}{4}} \\
s + \frac{1-s}{4} &= \frac{6}{23}\\
4s + (1-s) &= \frac{24}{23}\\
3s + 1 &=  \frac{24}{23} \\
3s &= \frac{24 - 23}{23} \\
s &= \frac{1}{23 \cdot 3} \approx 0.01449
\end{align}

Therefore approximately 1.45% of the code is serial.

