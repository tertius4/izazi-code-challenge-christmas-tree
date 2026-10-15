# Approach: Packing Christmas Trees

I wrote down my thought process, along with some remarks, throughout the project.

## What Am I Solving?
I started by thinking broadly: what do I want my solution to cover?
Clients regularly give feedback that changes the requirements. If your code is too specific, you may have to rewrite a large part of it. Sometimes two clients also want similar, but not identical, solutions. If your code is more general, you can reuse the same program for both.

I therefore wanted to support any polygon. The danger is that a too general approach might rule out some significant optimisations.
The instructions hinted that I should use y-symmetry. That helped me decide to write the program for symmetric polygons only.

## Which Algorithm?
While chatting back and forth with AI, it suggested using different algorithms for different numbers of trees (e.g. algorithm A for 1–60 trees and algorithm B for 61–100).
I liked the idea. It means more computation, but I thought that was worth it as long as the program finishes in a reasonable time.

## Which Programming Language?
C++ has the fastest computation.

I considered:
- **C** (lower level than C++): A quick Google search says it has no speed advantage. Fewer people know the language, and the code needs to be maintainable, so I did not opt for C.
- **Python (and JS):** These are 50–100 times slower than C++. Easy choice.
- **Fortran:** AI suggested it, but it seems better suited to pure mathematical computation.

The only catch is that I have to learn C++. I have experience in C, so it shouldn't be too difficult.

## Input & Output (I/O)
The output must be a .csv file (per the instructions). Since I wanted to support any symmetric polygon, I didn't want to hardcode the tree shape. Instead, the shape lives in a .csv file that the program reads.

## File Structure
My codebase was growing and I needed to organise it.
First, I made a GitHub repo and added the README.
Then I asked AI about standard file structures. Apparently there is no single standard for C++, but it suggested:
- `/src`: the program code
- `/build`: the compiled output
- `/data`: the input and output .csv files
- `/tests`: test cases

I figured the CLI could come later.

## The Body
I set the AI to "plan" mode and went back and forth on what a good solution would look like. During this dialogue, I realised that this problem came from a coding competition in 2025… lol.

Once the code was implemented, I wanted to be sure that I understood it (and that anybody else could too). I spent about an hour reading through it and changing small things. Someone is going to look at all my dirty laundry, and I want it clean!

## Evaluation
The visualisation is really cool. I liked that.
But my solution failed, due to overlap 😲. The visualisation showed no overlapping trees. Confusing...
So, I increased the spacing between them by a small fraction. After that, it passed :).

## Finishing Off
I added the CLI and ran the program a few times. I was not happy with the total score (it was ±95 at the time). AI suggested some changes and implemented them, which brought the score down to ±75. Good enough.
The downside is that the execution time on my machine is now up to 5 minutes :(
That sucks, but I don't think it's a big issue in production. Just go make a cup of coffee in the meantime.

This took me about 7 hours in total. To round off a "day's work", I went through the code once more to make sure it is easy to understand.

## Optimisation
I took a break, and thought about the code a bit. I think I called "Finishing off" too quickly. Surely, I should be able to do better than 75 in 5 mins...
