# Approach: Packing Christmas Trees

I wrote down my thought process, along with some remarks, throughout the project.

## Day 1

### What Am I Solving?
I started by thinking broadly: what do I want my solution to cover?
Clients change requirements, and two clients often want similar but not identical things. Choosing a more general approach works better in the long run; overly specific code usually has to be rewritten.

I therefore wanted to support any polygon. The danger is that too general an approach might rule out some significant optimisations.
The instructions hinted that I should use y-symmetry. That helped me decide to write the program for symmetric polygons only.

### Which Algorithm?
While chatting back and forth with AI, I got the suggestion to use different algorithms for different numbers of trees (e.g. algorithm A for 1–60 trees and algorithm B for 61–100).
I liked the idea. It means more computation, but I thought that was worth it as long as the program finishes in a reasonable time. I went for random starts for small numbers of trees and lattice patterns for large ones, then annealing for both.

### Which Programming Language?

I considered:
- **C++**: Seemed like a good option. Many people know it, so other programmers can maintain the code.
- **C**: A quick Google search says it has no speed advantage over C++. Fewer people know the language, so I did not opt for C.
- **Python (and JS):** These are 50–100 times slower than C++. Easy choice.
- **Fortran:** AI suggested it, but it seems better suited to pure mathematical computation.

The only catch is that I have to learn C++. I have experience in C, so it shouldn't be too difficult.

### Input & Output (I/O)
The output must be a .csv file (per the instructions). Since I wanted to support any symmetric polygon, I didn't want to hardcode the tree shape. Instead, the shape lives in a .csv file that the program reads.

### File Structure
My codebase was growing and I needed to organise it.
First, I made a GitHub repo and added the README.
Then I asked AI about standard file structures. Apparently there is no single standard for C++, but it suggested:
- `/src`: the program code
- `/build`: the compiled output
- `/data`: the input and output .csv files
- `/tests`: test cases

I figured a CLI could come later.

### The Algorithm
I set the AI to "plan" mode and went back and forth on what a good solution would look like. During this dialogue, I realised that this challenge came from a coding competition in 2025… lol.

Once the code was implemented, I wanted to be sure that I understood it (and that anybody else could too). I spent about an hour reading through it and changing small things. Someone is going to look at all my dirty laundry, and I want it clean!

How the final algorithm works is described in the [README](README.md#how-it-works).

### Evaluation
The visualisation is really cool. I liked that.
But my solution failed due to overlap 😲. The visualisation showed no overlapping trees. Confusing...
Turns out the trees overlapped by a rounding error far too small to see. So I added a tiny gap between them, and it passed 😃

### Finishing Off
I added the CLI and ran the program a few times. I was not happy with the total score (it was ±95 at the time). AI suggested some changes and implemented them, which brought the score down to ±75. Good enough.
The downside is that the execution time on my machine is now up to 5 minutes 😞
That sucks, but I don't think it's a big issue in production. Just go make a cup of coffee in the meantime.

This took me about 7 hours in total. I had to rush off to Bible study, so I decided to submit my result in the morning.

That night, while talking to friends, I thought: OK, I solved the problem, but is that the absolute best I can do? I have an entire week, after all.

## Day 2

### Quick Optimisation
In the morning I was eager to get the best score I could.
In quick succession, AI sped up packing 1–200 trees, then brought the score from ±75 → 73.75 → 72.

Curious where my result would place on Kaggle, I looked at the leaderboard. With 72, I'd be in ±1750th place.

<p align="center">
  <img src="./assets/disgusting.webp" alt="Disgusted reaction to my leaderboard position" width="300">
</p>

### No Modesty

It became clear that I had been thinking about this project wrong.

The goal was not to ship a good product to a client, but to get the absolute best score. I let go of all modesty – ignore run time, ignore the broad symmetric-polygon reasoning.

This is me and the tree! (…or the Christmas forest)

Goodbye, `tree.csv`: the tree is now hard-coded as four convex parts, which makes collision checks much faster.

### Trial and Error

This meant letting AI try approaches, testing them, discarding what did not work and keeping what did. So I started branching in GitHub, and had AI write a benchmark script for a sneak peek at whether a new strategy works, without waiting for all 200 packings, which by now took 50 minutes!

### Breakthrough

Shortly thereafter, it suggested keeping the best result in a .csv file and improving on it with every run: if a new result is better, it replaces the best; otherwise it's discarded.

This was what worked in the end. Now I am just constantly running the script on my computer to improve my result.

I stopped the script at: 70.xx.