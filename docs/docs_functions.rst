..
   ############################################################################
   ############################################################################
   ##                                                                        ##
   ##      ___  _   _   __   ____  ____   __   _   _ _____                   ##
   ##     / _ \| | | | /  \ |  _ \| __ \ /  \ | \ | |_   _|                  ##
   ##    ( |_| ) |_| |/ /\ \| |_| | -/ // /\ \|  \| | | |                    ##
   ##     \_  /|_____| /--\ |____/|_|\_\ /--\ |_\___| |_|                    ##
   ##       \/                                               v 0.0           ##
   ##                                                                        ##
   ############################################################################
   ############################################################################

.. image:: /_static/quadrant_logo.png

|

Function Reference 
==================

The order of functions in this API reference goes according to the chronological order of which they are called in the native LEOGPS processing work flow (see previous page). Functions that are currently not in use in the current native work flow are listed at the end of this page.

####

spacecraft.py
-------------

Spacecraft object

.. automodule:: spacecraft
   :members:

####

attitudes.py
------------

Attitude file containing attitude classes (quaternions, classical rodrigues parameters and modified rodrigues parameters).

.. automodule:: attitudes
   :members: QTR, CRP, MRP

####

This API reference was automatically generated using Sphinx' Autodoc feature, using the `NumPy docstring format <https://numpydoc.readthedocs.io/en/latest/format.html>`_, and last updated on 11th September 2021.
