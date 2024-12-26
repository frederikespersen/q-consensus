# Q-consensus: Read Quality Scores for Consensus Sequences

How do you interpret read Q-scores when you have several reads for the same sequence?
``Q-consensus`` is introduced to tackle that question.
``Q-consensus`` derives from the concept that more information is gained about a template sequence, the more reads that are provided; both when reads agree and disagree on a called base.

For more information, see ``draft_v0.3.pdf`` for a motivation for and derivation of the model.

**Note**: *Q-consensus* is a work in progress. See "To-do" section below.

## Contents

* ``draft_v0.3.pdf``: The derivation of the ``Q-consensus`` model [Draft]
* ``q-consensus.py``: A Python implementation of the ``Q-consensus`` model
* ``read-depth.ipynb``: Testing how the model predicts Q-scores to increase when multiple reads agree on a consensus.

## To-do

* Handle gaps in alignment (and final refinement of model)
* Optimization of speed
* Set up as python module
* Third party validation
