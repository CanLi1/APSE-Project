## Duran 8 Process Problem
This is a re-implementation of Duran example 3 superstructure synthesis (the "eight process problem") in Pyomo with outer approximation cuts generated automatically by automatic differentiation.

Ref:
    SELECT OPTIMAL PROCESS FROM WITHIN GIVEN SUPERSTRUCTURE.
    MARCO DURAN , PH.D. THESIS (EX3) , 1984.
    CARNEGIE-MELLON UNIVERSITY , PITTSBURGH , PA.

(original problem, my implementation may vary)
    Problem type:   convex MINLP
            size:    8  binary variables
                    26  continuous variables
                    32  constraints

Pictoral representation can be found on
page 969 of [Turkay & Grossmann, 1996](http://dx.doi.org/10.1016/0098-1354(95)00219-7).

Used result of set covering algorithm from Turkay paper cited above for initial NLP selection

Two variations are currently provided:
* main.py - This is the classic variation with a single-component system
* main_2comp.py - This is an extension for two components, meant primarily as a test-bed for handling of indexed variables in the automatic differentiation process. As the splitters do not currently enforce bilinear split fraction constraints, this is actually a relaxation of the original problem.