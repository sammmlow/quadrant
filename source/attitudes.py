# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:18:36 2021

@author: sammm
"""

import numpy as np

# Class for quarternions
class QRT:
    
    ########################################################################
    ########################################################################
    
    def __init__(self, q=None, dcm=None):
        
        # If the coordinates and the DCM are not specified, then just return
        # an quarternion and DCM with zero rotation.
        if q is None and dcm is None:
            self.q = np.array([1.0, 0.0, 0.0, 0.0])
            self.dcm = np.eye(3)
        
        # If the quarternion is specified, regardless of what the specified
        # DCM is, the quarternion will take precedence.
        elif q is not None:
            if len(q) == 4:
                self.q = np.array(q)
                self.dcm = self._qrt2dcm( self.q )
            else:
                print("Error! Quarternion not of length 4!")
                self.q = np.array([1.0, 0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
        
        # If the quarternion is not specified, but the DCM is, convert the DCM
        # into a quarternion and return both attributes.
        elif q is None and dcm is not None:
            if np.shape(dcm)[0] == 3 and np.shape(dcm)[1] == 3:
                self.q = self._dcm2qrt( dcm )
                self.dcm = dcm
            else:
                print("Error! Quarternion not specified and DCM is not 3x3!")
                self.q = np.array([1.0, 0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
    
    ########################################################################
    ########################################################################
    
    # Property setting.
    
    @property
    def q(self):
        return self._q
    
    @q.setter
    def q(self, value):
        self._q = value
        self.dcm = self._qrt2dcm( value )
    
    ########################################################################
    ########################################################################
    
    def __add__(self, summand):
        
        # If there was a rotation FN, done by sequentially [FB]*[BN]
        # Then the summand represents quarternions in the [BN] DCM
        # And the self represents quarternions in the [FB] DCM
        # Note that quarternion addition is not commutative!
        
        qC  = self.q
        qC2 = summand.q
        
        qM = np.array([[qC[0], -1*qC[1], -1*qC[2], -1*qC[3]],
                       [qC[1],    qC[0],    qC[3], -1*qC[2]],
                       [qC[2], -1*qC[3],    qC[0],    qC[1]],
                       [qC[3],    qC[2], -1*qC[1],    qC[0]]])
        
        # Now we can return the quarternion [FB].
        return QRT( q=np.transpose( qM @ np.transpose(qC2) ) )
    
    ########################################################################
    ########################################################################
    
    def __sub__(self, minuend):
        
        # If there was a rotation DCM [FN], done by sequentially [FB]*[BN]
        # Then the minuend represents quarternions in the [BN] DCM
        # And the self represents the final quarternion in the [FN] DCM
        # Thus, if the user specifies qOut = q2 - q1, then the inputs follow
        # the order such that [FN( q2 )] = [FB( qOut )] * [BN( q1 )].
        # Note that quarternion subtraction is not commutative!
        
        qC, qC2  = self.q, np.zeros(4)
        
        # Flip the angular components of the minuend quarternion.
        qC2[0] =    minuend.q[0]
        qC2[1] = -1*minuend.q[1]
        qC2[2] = -1*minuend.q[2]
        qC2[3] = -1*minuend.q[3]
        
        qM = np.array([[qC[0], -1*qC[1], -1*qC[2], -1*qC[3]],
                       [qC[1],    qC[0],    qC[3], -1*qC[2]],
                       [qC[2], -1*qC[3],    qC[0],    qC[1]],
                       [qC[3],    qC[2], -1*qC[1],    qC[0]]])
        
        # Now we have [FB] = [FN]*[NB]
        return QRT( q = np.transpose( qM @ np.transpose(qC2) ) )
    
    ########################################################################
    ########################################################################
    
    def __iter__(self):
        return iter(self.q)
    
    ########################################################################
    ########################################################################
    
    def __next__(self):
        return next(self.q)
    
    ########################################################################
    ########################################################################
    
    def __getitem__(self, obj):
        return self.q[obj]
    
    ########################################################################
    ########################################################################
    
    # Representation
    def __repr__(self):
        return repr(self.q)
    
    ########################################################################
    ########################################################################
    
    # Representation if returned as a string
    def __str__(self):
        return str(self.q)
    
    ########################################################################
    ########################################################################
    
    # Length of quarternion vector.
    def __len__(self):
        return len(self.q)
    
    ########################################################################
    ########################################################################
    
    # Method for conversion of a direction cosine matrix into a quarternion.
    def _dcm2qrt(self, dcm):
    
        # First, calculate the individual betas of the quarternion vector.
        b0sq = 0.25 * ( 1 + np.trace(dcm) )
        b1sq = 0.25 * ( 1 + 2*dcm[0,0] - np.trace(dcm) )
        b2sq = 0.25 * ( 1 + 2*dcm[1,1] - np.trace(dcm) )
        b3sq = 0.25 * ( 1 + 2*dcm[2,2] - np.trace(dcm) )
        
        # Next, calculate the intermediate betas used in Sheppard's method.
        b0b1 = ( dcm[1,2] - dcm[2,1] ) / 4
        b0b2 = ( dcm[2,0] - dcm[0,2] ) / 4
        b0b3 = ( dcm[0,1] - dcm[1,0] ) / 4
        b1b2 = ( dcm[0,1] + dcm[1,0] ) / 4
        b3b1 = ( dcm[2,0] + dcm[0,2] ) / 4
        b2b3 = ( dcm[1,2] + dcm[2,1] ) / 4
        
        # Check to see which beta is largest.
        bmax = max( [ b0sq, b1sq, b2sq, b3sq ] )
        
        # Divide by the largest beta.
        if bmax == b0sq:
            b0 = b0sq ** 0.5
            b1 = b0b1 / b0
            b2 = b0b2 / b0
            b3 = b0b3 / b0
            
        elif bmax == b1sq:
            b1 = b1sq ** 0.5
            b0 = b0b1 / b1
            b2 = b1b2 / b1
            b3 = b3b1 / b1
            
        elif bmax == b2sq:
            b2 = b2sq ** 0.5
            b0 = b0b2 / b2
            b1 = b1b2 / b2
            b3 = b2b3 / b2
            
        elif bmax == b3sq:
            b3 = b3sq ** 0.5
            b0 = b0b3 / b3
            b1 = b3b1 / b3
            b2 = b2b3 / b3
        
        # Output a 1x4 Quarternion
        self.q = np.array([b0, b1, b2, b3])
        return self.q
    
    ########################################################################
    ########################################################################
    
    # Method for conversion of a quarternion into a direction cosine matrix.
    def _qrt2dcm(self, qrt):
        
        # First, get the individual quarternion components.
        b0, b1, b2, b3 = qrt[0], qrt[1], qrt[2], qrt[3]
        
        # Next, get the direction cosine matrix.
        dcm = np.zeros((3,3))
        
        # Compute the DCM diagonals
        dcm[0,0] = b0**2 + b1**2 - b2**2 - b3**2
        dcm[1,1] = b0**2 - b1**2 + b2**2 - b3**2
        dcm[2,2] = b0**2 - b1**2 - b2**2 + b3**2
        
        # Compute the DCM off-diagonals
        dcm[0,1] = 2 * (b1 * b2 + b0 * b3)
        dcm[0,2] = 2 * (b1 * b3 - b0 * b2)
        dcm[1,0] = 2 * (b1 * b2 - b0 * b3)
        dcm[1,2] = 2 * (b2 * b3 + b0 * b1)
        dcm[2,0] = 2 * (b1 * b3 + b0 * b2)
        dcm[2,1] = 2 * (b2 * b3 - b0 * b1)
        
        self.dcm = dcm
        return self.dcm
    
    ########################################################################
    ########################################################################
    
    # Method for computing the differential equation in quarternion vector
    # form given a quarternion and an angular velocity vector.
    def compute_qrate(self, omega):
        
        # Input: 1x4 quarternion coordinates and 1x3 angular velocity vector.
        # Returns the 1x4 quarternion dot (derivative) vector.
        
        # Note that the transformation matrix qM that relates d(beta)/dt and
        # the angular velocity vector is orthogonal and singularity free.
        # This means the inverse transformation can also be defined.
        
        qM = np.array([[-1*self.q[1], -1*self.q[2], -1*self.q[3]],
                       [   self.q[0], -1*self.q[3],    self.q[2]],
                       [   self.q[3],    self.q[0], -1*self.q[1]],
                       [-1*self.q[2],    self.q[1],    self.q[0]]])
        
        self.qrate = 0.5 * np.transpose( qM @ omega )
        return self.qrate