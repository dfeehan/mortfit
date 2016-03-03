**mortfit**
-------

Fit mortality hazard models, with special attention to humans at old ages.

Package currently under development. Stay tuned.

TODO
----

- go through and remove subset(); this should only be used interactively. see https://github.com/hadley/devtools/wiki/Evaluation for useful info


ideas for tests
---------------

* be sure that the predicted values from all of the pre-defined functional forms always make sense (probs between 0 and 1, numbers of deaths never more than number of people, etc)

* be sure that optim always converges on test datasets

* test for haz.to.prob: one way of doing this is to take fairly easy cases, with known closed-form solutions, and to test the haz.to.prob output against them



