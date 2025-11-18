---
title: "Assignment 04"
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

# Assignment 4 (Poison Solver)

## (a) Scalability Study on One Node

### Domain Size: 100×100 (omega = 1.1)

| Ranks | STD Time (s) | RB Time (s) | STD Speedup | RB Speedup | RB Advantage |
|-------|--------------|-------------|-------------|------------|--------------|
| 1     | 2.82         | 0.59        | 1.00×       | 1.00×      | 4.78×        |
| 2     | 1.47         | 0.36        | 1.92×       | 1.64×      | 4.08×        |
| 4     | 0.79         | 0.25        | 3.57×       | 2.36×      | 3.16×        |
| 8     | 0.47         | 0.20        | 6.00×       | 2.95×      | 2.35×        |
| 16    | 0.30         | 0.16        | 9.40×       | 3.69×      | 1.88×        |
| 32    | 0.25         | 0.17        | 11.28×      | 3.47×      | 1.47×        |
| 64    | 0.35         | 0.30        | 8.06×       | 1.97×      | 1.17×        |
| 72    | 0.45         | 0.43        | 6.27×       | 1.37×      | 1.05×        |


**Observations**: Red-Black  is consistently faster. Both solvers show best performance around 16-32 processes Afterwards degradation beyond 32 ranks due to  communication overhead dominates for small domains

### Domain Size: 400×400 (omega = 1.1)

| Ranks | STD Time (s) | RB Time (s) | STD Speedup | RB Speedup | RB Advantage |
|-------|--------------|-------------|-------------|------------|--------------|
| 1     | 756.11       | 148.72      | 1.00×       | 1.00×      | 5.08×        |
| 2     | 378.94       | 75.81       | 2.00×       | 1.96×      | 5.00×        |
| 4     | 190.72       | 39.20       | 3.96×       | 3.79×      | 4.87×        |
| 8     | 96.66        | 20.97       | 7.82×       | 7.09×      | 4.61×        |
| 16    | 49.81        | 11.93       | 15.18×      | 12.47×     | 4.18×        |
| 32    | 27.21        | 7.45        | 27.79×      | 19.96×     | 3.65×        |
| 64    | 16.79        | 5.94        | 45.03×      | 25.04×     | 2.83×        |
| 72    | 15.21        | 5.58        | 49.71×      | 26.66×     | 2.73×        |



**Observations**: Near-linear scaling up to 32 ranks. RB maintains is consistent approx. 3-5 times faster. Communication overhead becomes less critical.

## (b) Solver Stability Analysis

### Impact of Domain Size (omega = 1.1)

| Domain | Ranks | STD Iterations | RB Iterations | Iteration Increase |
|--------|-------|----------------|---------------|-------------------|
| 100    | 1     | 26,867         | 26,870        | ~0%               |
| 200    | 1     | 106,635        | 106,639       | ~0%               |
| 300    | 1     | 239,289        | 239,293       | ~0%               |
| 400    | 1     | 424,829        | 424,833       | ~0%               |

**Observation**:Iterations scale approximately as O of N squared, where N is the domain size.

### Impact of Process Count (Domain = 250×250, omega = 1.1)

| Ranks | STD Iterations | RB Iterations | STD Degradation | RB Degradation |
|-------|----------------|---------------|-----------------|----------------|
| 1     | 166,352        | 166,355       | baseline        | baseline       |
| 2     | 166,586        | 166,355       | +0.14%          | 0%             |
| 4     | 167,268        | 166,302       | +0.55%          | -0.03%         |
| 8     | 168,873        | 166,383       | +1.52%          | +0.02%         |
| 16    | 171,278        | 165,721       | +2.96%          | -0.38%         |
| 32    | 177,735        | 166,059       | +6.84%          | -0.18%         |
| 64    | 189,025        | 165,400       | +13.62%         | -0.57%         |
| 72    | 190,416        | 164,079       | +14.46%         | -1.37%         |

**Observation**: STD solver degrades significantly and has to use up to 14.5% more iterations needed at 72 ranks and RB solver remains stable: Iteration count slightly reduced with higher number of ranks

### Impact of Omega (Domain = 250×250)

#### omega = 1.5 (Moderate Over-relaxation)

| Ranks | STD Iterations | RB Iterations | STD Convergence | RB Convergence |
|-------|----------------|---------------|-----------------|----------------|
| 1     | 70,937         | 70,940        | Yes             | Yes            |
| 16    | 76,454         | 70,688        | Yes             | Yes            |
| 32    | 82,928         | 70,836        | Yes             | Yes            |
| 64    | 94,959         | 70,553        | Yes             | Yes            |
| 72    | 97,266         | 70,014        | Yes             | YES            |

**Observation**: Both solvers are stable, but STD shows stronge degradation in comparison to RB.

#### omega = 1.8 (Aggressive Over-relaxation)

| Ranks | STD Iterations | RB Iterations | STD Convergence | RB Convergence |
|-------|----------------|---------------|-----------------|----------------|
| 1     | 25,741         | 25,741        | Yes             | Yes            |
| 16    | 31,743         | 25,665        | Yes             | Yes            |
| 32    | 20,916         | 25,726        | **NaN**       | Yes            |
| 64    | 3,437          | 25,615        | **NaN**       | Yes            |
| 72    | 3,039          | 25,435        | **NaN**       | Yes            |

**Observation**: STD solver becomes unstable at high rank counts and with aggressive over-relaxation.

#### ω = 1.95 (Near-instability threshold)

| Ranks | STD Iterations | RB Iterations | STD Convergence | RB Convergence |
|-------|----------------|---------------|-----------------|----------------|
| 1     | 7,340          | 7,331         | Yes             | Yes            |
| 4     | 8,585          | 7,333         | Yes             | Yes            |
| 8     | 18,781         | 7,335         | **NaN**       | Yes            |
| 16    | 7,487          | 7,324         | **NaN**       | Yes            |
| 32    | 3,756          | 7,355         | **NaN**       | Yes            |
| 64    | 1,961          | 7,308         | **NaN**       | Yes            |
| 72    | 1,805          | 7,270         | **NaN**       | Yes            |

**Observation**: STD becomes even more unstable at earlier ranks.

## (c) Convergence Comparison to Sequential Case

### Iteration Count Analysis

| Domain | omega    | Sequential STD | Sequential RB | Parallel STD (64 ranks) | Parallel RB (64 ranks) |
|--------|------|----------------|---------------|-------------------------|------------------------|
| 100    | 1.1  | 26,867         | 26,870        | 35,918 (+33.7%)         | 26,458 (-1.5%)         |
| 200    | 1.1  | 106,635        | 106,639       | 125,675 (+17.9%)        | 107,086 (+0.4%)        |
| 300    | 1.1  | 239,289        | 239,293       | 263,647 (+10.2%)        | 235,603 (-1.5%)        |
| 400    | 1.1  | 424,829        | 424,833       | 461,916 (+8.7%)         | 424,516 (-0.1%)        |
| 250    | 1.5  | 70,937         | 70,940        | 94,959 (+33.9%)         | 70,553 (-0.5%)         |
| 250    | 1.8  | 25,741         | 25,741        | **NaN (diverged)**      | 25,615 (-0.5%)         |

**Observations**: RB maintains number of iterations in sequential convergence. STD degrades significantly with 8-34% more iterations needed in parallel. Degradation is worse for small domains.

## (d) Response to Different Omega Values

### Optimal Omega Analysis

| omega    | Description           | STD Sequential | RB Sequential | STD Parallel (64) | RB Parallel (64) |
|------|-----------------------|----------------|---------------|-------------------|------------------|
| 1.1  | Safe under-relaxation | 166,352 iter   | 166,355 iter  | 189,025 iter      | 165,400 iter     |
| 1.5  | Moderate SOR          | 70,937 iter    | 70,940 iter   | 94,959 iter       | 70,553 iter      |
| 1.8  | Aggressive SOR        | 25,741 iter    | 25,741 iter   | **Diverges**      | 25,615 iter      |
| 1.95 | Near-optimal          | 7,340 iter     | 7,331 iter    | **Diverges**      | 7,308 iter       |

The **Optimal omega** is around 1.95 for this problem.
