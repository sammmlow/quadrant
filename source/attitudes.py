# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:18:36 2021

@author: sammm
"""

import numpy as np


###############################################################################
###############################################################################
##                                                                           ##
##                 Attitude Coordinate: Quaternion Class (QTR)               ##
##                                                                           ##
###############################################################################
###############################################################################


class QTR:
    '''
    A Quaternion (QTR) class can be initialised with either a 1x4
    NumPy array, or a direction cosine as a 3x3 NumPy array.
    
    Parameters
    ----------
    qtr : TYPE, optional
        DESCRIPTION. The default is None.
    dcm : TYPE, optional
        DESCRIPTION. The default is None.

    '''
    
    def __init__(self, qtr=None, dcm=None):
        
        
        # If the coordinates and the DCM are not specified, then just return
        # a quaternion and DCM with zero rotation.
        if qtr is None and dcm is None:
            self.qtr = np.array([1.0, 0.0, 0.0, 0.0])
            self.dcm = np.eye(3)
        
        # If the quaternion is specified, regardless of what the specified
        # DCM is, the quaternion will take precedence.
        elif qtr is not None:
            if len(qtr) == 4:
                self.qtr = np.array(qtr)
                self.normalise()
                self.dcm = self._qtr2dcm( self.qtr )
            else:
                print("Error! Quaternion not of length 4!")
                raise ValueError("Quaternion not of length 4!")
                self.qtr = np.array([1.0, 0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
        
        # If the quaternion is not specified, but the DCM is, convert the DCM
        # into a quaternion and return both attributes.
        elif qtr is None and dcm is not None:
            if np.shape(dcm)[0] == 3 and np.shape(dcm)[1] == 3:
                self.qtr = self._dcm2qtr( dcm )
                self.normalise()
                self.dcm = dcm
            else:
                print("Error! QTR not specified and DCM is not 3x3!")
                self.qtr = np.array([1.0, 0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
    
    # Property setting so that any updates to self.q updates self.dcm
    @property
    def qtr(self):
        return self._qtr
    @qtr.setter
    def qtr(self, value):
        self._qtr = value
        self._qtr2dcm( value )
    
    # Quaternion product. Note that quaternion product is not commutative!
    # Self should be the parent rotation (qP), and summand is the child (qC).
    # qF = qP * qC ---> qP is applied first, then qC second.
    def __mul__(self, summand):
        qP, qC = self.qtr, summand.qtr
        qM = np.array([[qC[0], -1*qC[1], -1*qC[2], -1*qC[3]],
                       [qC[1],    qC[0],    qC[3], -1*qC[2]],
                       [qC[2], -1*qC[3],    qC[0],    qC[1]],
                       [qC[3],    qC[2], -1*qC[1],    qC[0]]])
        return QTR( qtr=np.transpose( qM @ qP ) )
    
    # Quaternion quotient. Note that quaternion quotient is not commutative!
    # Self should be the parent rotation (qP), and minuend is the child (qC).
    # qF = qP / qC ---> qP is applied first, then qC inverse is second.
    def __truediv__(self, minuend):
        qP = self.qtr
        qC = np.append( minuend[0], -1*minuend[1:] ) 
        qM = np.array([[qC[0], -1*qC[1], -1*qC[2], -1*qC[3]],
                       [qC[1],    qC[0],    qC[3], -1*qC[2]],
                       [qC[2], -1*qC[3],    qC[0],    qC[1]],
                       [qC[3],    qC[2], -1*qC[1],    qC[0]]])
        return QTR( qtr = np.transpose( qM @ qP ) )
    
    # String ID.
    def strID(self):
        return 'QTR'
    
    # Iterator.
    def __iter__(self):
        return iter(self.qtr)
    
    # Next attribute for iteration.
    def __next__(self):
        return next(self.qtr)
    
    # Get item.
    def __getitem__(self, obj):
        return self.qtr[obj]
    
    # String representation.
    def __repr__(self):
        return "attitude_qtr" + repr(list(self.qtr))
    
    # Representation if returned as a string
    def __str__(self):
        return str(self.qtr)
    
    # Length of quaternion vector.
    def __len__(self):
        return len(self.qtr)
    
    # Method for conversion of a direction cosine matrix into a quaternion.
    def _dcm2qtr(self, dcm):
    
        # First, calculate the individual betas of the quaternion vector.
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
        
        self.qtr = np.array([b0, b1, b2, b3])
        return self.qtr
    
    # Method for conversion of a quaternion into a direction cosine matrix.
    def _qtr2dcm(self, qtr):
        
        # First, get the individual quaternion components.
        b0, b1, b2, b3 = qtr[0], qtr[1], qtr[2], qtr[3]
        
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
    
    # Method for ensuring the correct normalisation of a quaternion.
    def normalise(self):
        self.qtr = self.qtr / np.linalg.norm(self.qtr)
    
    # Method for computing the differential equation in quaternion vector
    # form given a quaternion and an angular velocity vector. Note that the
    # attitude rate is NOT a property of the attitude coordinate instance as
    # it relies on the angular velocity of the spacecraft, which is external.
    def get_qtrRate(self, omega):
        
        # Note that the transformation matrix qM that relates d(beta)/dt and
        # the angular velocity vector is orthogonal and singularity free.
        # This means the inverse transformation can also be defined.
        qM = np.array([[-1*self.qtr[1], -1*self.qtr[2], -1*self.qtr[3]],
                       [   self.qtr[0], -1*self.qtr[3],    self.qtr[2]],
                       [   self.qtr[3],    self.qtr[0], -1*self.qtr[1]],
                       [-1*self.qtr[2],    self.qtr[1],    self.qtr[0]]])
        
        rate = 0.5 * np.transpose( qM @ omega )
        return rate
    
    # Method for obtaining the angle of rotation in degrees described by
    # this attitude coordinate.
    def get_angle(self):
        return 2 * np.rad2deg( np.arccos( self.qtr[0] ) )
    
    # Method for obtaining the axis of rotation that this attitude describes.
    # Note that if the attitude describes a transformation from the B frame 
    # to the N frame, or vice versa, the axis of rotation is the same in both
    # frames and which is why it is an eigenvector of the rotation matrix.
    def get_axis(self):
        axis = np.array([ self.qtr[1], self.qtr[2], self.qtr[3] ])
        axis = axis / np.linalg.norm(axis)
        return axis


###############################################################################
###############################################################################
##                                                                           ##
##       Attitude Coordinate: Classical Rodrigues Parameter Class (CRP)      ##
##                                                                           ##
###############################################################################
###############################################################################


class CRP:
    
    def __init__(self, crp=None, dcm=None):
        '''
        A Classical Rodrigues Parameter class can be initialised with 
        either a 1x3 NumPy array, or a direction cosine as a 3x3 NumPy array.
        
        Parameters
        ----------
        crp : TYPE, optional
            DESCRIPTION. The default is None.
        dcm : TYPE, optional
            DESCRIPTION. The default is None.

        '''
        
        # If the coordinates and the DCM are not specified, then just return
        # a CRP and DCM with zero rotation.
        if crp is None and dcm is None:
            self.crp = np.array([0.0, 0.0, 0.0])
            self.dcm = np.eye(3)
        
        # If the CRP is specified, regardless of what the specified
        # DCM is, the CRP will take precedence.
        elif crp is not None:
            if len(crp) == 3:
                self.crp = np.array(crp)
                self.dcm = self._crp2dcm( self.crp )
            else:
                print("Error! Classical Rodrigues Parameter not of length 3!")
                raise ValueError("CRP not of length 3!")
                self.crp = np.array([0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
        
        # If the CRP is not specified, but the DCM is, convert the DCM
        # into a CRP and return both attributes.
        elif crp is None and dcm is not None:
            if np.shape(dcm)[0] == 3 and np.shape(dcm)[1] == 3:
                self.crp = self._dcm2crp( dcm )
                self.dcm = dcm
            else:
                print("Error! CRP not specified and DCM is not 3x3!")
                self.crp = np.array([0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
    
    # Property setting so that any updates to self.q updates self.dcm
    @property
    def crp(self):
        return self._crp
    @crp.setter
    def crp(self, value):
        self._crp = value
        self._crp2dcm( value )
    
    # CRP product. Note that the CRP product is not commutative!
    # Self should be the parent rotation (cP), and summand is the child (cC).
    # cF = cP * cC ---> cP is applied first, then cC second.
    def __mul__(self, summand):
        cC, cP = self.crp, summand.crp
        return CRP(crp=((cP + cC - (np.cross(cP,cC))) / (1 - np.dot(cP,cC))))
    
    # CRP quotient. Note that the CRP quotient is not commutative!
    # Self should be the parent rotation (cP), and minuend is the child (cC).
    # cF = cP / cC ---> cP is applied first, then cC inverse is second.  
    def __truediv__(self, minuend):
        cC, cP = self.crp, minuend.crp
        return CRP(crp=((cC - cP + np.cross(cP,cC)) / (1 + np.dot(cP,cC))))
    
    # String ID.
    def strID(self):
        return 'CRP'
    
    # Iterator.
    def __iter__(self):
        return iter(self.crp)
    
    # Next attribute for iteration.
    def __next__(self):
        return next(self.crp)
    
    # Get item.
    def __getitem__(self, obj):
        return self.crp[obj]
    
    # String representation.
    def __repr__(self):
        return "attitude_crp" + repr(list(self.crp))
    
    # Representation if returned as a string
    def __str__(self):
        return str(self.crp)
    
    # Length of CRP vector.
    def __len__(self):
        return len(self.crp)
    
    # Method for conversion of a direction cosine matrix into a CRP.
    def _dcm2crp(self, dcm):
        c1 = dcm[1,2] - dcm[2,1]
        c2 = dcm[2,0] - dcm[0,2]
        c3 = dcm[0,1] - dcm[1,0]
        self.crp = ( 1 / ( np.trace(dcm) + 1 ) ) * np.array([c1,c2,c3])
        return self.crp
    
    # Method for conversion of a CRP into a direction cosine matrix.
    def _crp2dcm(self, crp):
        
        # Initialise the DCM and its diagonals.
        dcm = np.zeros((3,3))
        dcm[0,0] = 1 + crp[0]**2 - crp[1]**2 - crp[2]**2
        dcm[1,1] = 1 - crp[0]**2 + crp[1]**2 - crp[2]**2
        dcm[2,2] = 1 - crp[0]**2 - crp[1]**2 + crp[2]**2
        
        # Initialise the DCM off-diagonals
        dcm[0,1] = 2 * (crp[0] * crp[1] + crp[2])
        dcm[0,2] = 2 * (crp[0] * crp[2] - crp[1])
        dcm[1,0] = 2 * (crp[1] * crp[0] - crp[2])
        dcm[1,2] = 2 * (crp[1] * crp[2] + crp[0])
        dcm[2,0] = 2 * (crp[2] * crp[0] + crp[1])
        dcm[2,1] = 2 * (crp[2] * crp[1] - crp[0])
        
        # Return the scaled DCM.
        self.dcm = ( 1.0 / ( 1.0 + np.dot(crp,crp) ) ) * dcm
        return self.dcm
        
    # Method for computing the differential equation in CRP vector form given
    # the self CRP and an angular velocity vector. Note that the attitude
    # rate is NOT a property of the attitude coordinate instance as it relies
    # on the angular velocity of the spacecraft, which is external.
    def get_crpRate(self, omega):
        
        # Compute the derivative transformation matrix.
        crp_matrix = np.zeros((3,3))
        crp_matrix[0,0] = 1 + self.crp[0]**2
        crp_matrix[0,1] = self.crp[0] * self.crp[1] - self.crp[2]
        crp_matrix[0,2] = self.crp[0] * self.crp[2] + self.crp[1]
        crp_matrix[1,0] = self.crp[1] * self.crp[0] + self.crp[2]
        crp_matrix[1,1] = 1 + self.crp[1]**2
        crp_matrix[1,2] = self.crp[1] * self.crp[2] - self.crp[0]
        crp_matrix[2,0] = self.crp[2] * self.crp[0] - self.crp[1]
        crp_matrix[2,1] = self.crp[2] * self.crp[1] + self.crp[0]
        crp_matrix[2,2] = 1 + self.crp[2]**2
        
        # Return the CRP derivative vector.
        rate = 0.5 * np.transpose( crp_matrix @ omega )
        return rate
    
    # Method for obtaining the angle of rotation in degrees described by
    # this attitude coordinate.
    def get_angle(self):
        return 2 * np.rad2deg( np.arctan( (np.linalg.norm( self.crp )) ) )
    
    # Method for obtaining the axis of rotation that this attitude describes.
    # Note that if the attitude describes a transformation from the B frame 
    # to the N frame, or vice versa, the axis of rotation is the same in both
    # frames and which is why it is an eigenvector of the rotation matrix.
    def get_axis(self):
        axis = np.array([ self.crp[0], self.crp[1], self.crp[2] ])
        axis = axis / np.linalg.norm(axis)
        return axis


###############################################################################
###############################################################################
##                                                                           ##
##       Attitude Coordinate: Modified Rodrigues Parameter Class (MRP)       ##
##                                                                           ##
###############################################################################
###############################################################################


class MRP:
    
    def __init__(self, mrp=None, dcm=None):
        '''
        A Modified Rodrigues Parameter class can be initialised with 
        either a 1x3 NumPy array, or a direction cosine as a 3x3 NumPy array.
        
        Parameters
        ----------
        mrp : TYPE, optional
            DESCRIPTION. The default is None.
        dcm : TYPE, optional
            DESCRIPTION. The default is None.

        '''
        
        # If the coordinates and the DCM are not specified, then just return
        # an MRP and DCM with zero rotation.
        if mrp is None and dcm is None:
            self.mrp = np.array([0.0, 0.0, 0.0])
            self.dcm = np.eye(3)
        
        # If the MRP is specified, regardless of what the specified
        # DCM is, the MRP will take precedence.
        elif mrp is not None:
            if len(mrp) == 3:
                self.mrp = np.array(mrp)
                self.dcm = self._mrp2dcm( self.mrp )
            else:
                print("Error! Modified Rodrigues Parameter not of length 3!")
                raise ValueError("MRP not of length 3!")
                self.mrp = np.array([0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
        
        # If the MRP is not specified, but the DCM is, convert the DCM
        # into a MRP and return both attributes.
        elif mrp is None and dcm is not None:
            if np.shape(dcm)[0] == 3 and np.shape(dcm)[1] == 3:
                self.mrp = self._dcm2mrp( dcm )
                self.dcm = dcm
            else:
                print("Error! MRP not specified and DCM is not 3x3!")
                self.mrp = np.array([0.0, 0.0, 0.0])
                self.dcm = np.eye(3)
    
    # Property setting so that any updates to the coordinate updates self.dcm.
    @property
    def mrp(self):
        return self._mrp
    @mrp.setter
    def mrp(self, value):
        self._mrp = value
        self._mrp2dcm( value )
    
    # MRP product. Note that the MRP product is not commutative!
    # Self should be the parent rotation (mP), and summand is the child (mC).
    # mF = mP * mC ---> mP is applied first, then mC second.
    def __mul__(self, summand):
        mC, mP = self.mrp, summand.mrp
        mCx = 1.0 - np.dot(mC,mC)
        mPx = 1.0 - np.dot(mP,mP)
        mrpi = (mCx * mP) + (mPx * mC) - 2 * np.cross(mP,mC)
        mrpi = mrpi / (1+(np.dot(mC,mC) * np.dot(mP,mP)) - (2*np.dot(mC,mP)))
        return MRP(mrp=mrpi)
    
    # MRP quotient. Note that the MRP quotient is not commutative!
    # Self should be the parent rotation (mP), and minuend is the child (mC).
    # mF = mP / mC ---> mP is applied first, then mC inverse is second.  
    def __truediv__(self, minuend):
        mC, mP = self.mrp, minuend.mrp
        mCx = 1.0 - np.dot(mC,mC)
        mPx = 1.0 - np.dot(mP,mP)
        mrpi = (mPx * mC) - (mCx * mP) + 2 * np.cross(mP,mC)
        mrpi = mrpi / (1+(np.dot(mC,mC) * np.dot(mP,mP)) + (2*np.dot(mC,mP)))
        return MRP(mrp=mrpi)
    
    # String ID.
    def strID(self):
        return 'MRP'
    
    # Iterator.
    def __iter__(self):
        return iter(self.mrp)
    
    # Next attribute for iteration.
    def __next__(self):
        return next(self.mrp)
    
    # Get item.
    def __getitem__(self, obj):
        return self.mrp[obj]
    
    # String representation.
    def __repr__(self):
        return "attitude_mrp" + repr(list(self.mrp))
    
    # Representation if returned as a string
    def __str__(self):
        return str(self.mrp)
    
    # Length of MRP vector.
    def __len__(self):
        return len(self.mrp)
    
    def _mrp_shadow(self):
        self.mrp = -1 * self.mrp / np.dot(self.mrp,self.mrp)
        return self.mrp
    
    # Method for conversion of a direction cosine matrix into a MRP. This
    # function will always provide the inner MRP set where |mrp| <= 1.
    def _dcm2mrp(self, dcm):
        zeta = (np.trace(dcm) + 1)**0.5
        
        # Check if a singularity occurs for 180 degree rotation.
        if zeta != 0.0:
            mt = (np.transpose(dcm)-dcm) / (zeta*(zeta+2))
            self.mrp = np.array([ mt[2,1], mt[0,2], mt[1,0]])
            
        # If there was a singularity, convert the DCM to quaternion.
        else:
            q = QTR( dcm=dcm )
            self.mrp = np.array([q[1], q[2], q[3]]) / (1 + q[0])
        
        return self.mrp
    
    # Method for conversion of a MRP into a direction cosine matrix.
    def _mrp2dcm(self, mrp):
        
        # Express the MRP with cross-product operator in tilde form.
        mt = np.array([[0.0, -1*mrp[2], mrp[1]],
                       [mrp[2], 0.0, -1*mrp[0]],
                       [-1*mrp[1], mrp[0], 0.0]])
        
        # Now we can form the DCM.
        dcm = 8 * (mt @ mt) - 4 * ( 1 - np.dot(mrp,mrp) ) * mt
        dcm = dcm / ((1 + np.dot(mrp,mrp))**2)
        dcm = dcm + np.eye(3)
        self.dcm = dcm
        return self.dcm
        
    # Method for computing the differential equation in MRP vector form given
    # the self MRP and an angular velocity vector. Note that the attitude
    # rate is NOT a property of the attitude coordinate instance as it relies
    # on the angular velocity of the spacecraft, which is external.
    def get_mrpRate(self, omega):
        
        # Express the MRP with cross-product operator in tilde form.
        mt = np.array([[0.0, -1*self.mrp[2], self.mrp[1]],
                       [self.mrp[2], 0.0, -1*self.mrp[0]],
                       [-1*self.mrp[1], self.mrp[0], 0.0]])
        
        # Compute the derivative transformation matrix B.
        B = (( 1 - np.dot(self.mrp, self.mrp) ) * np.eye(3))
        B = B + (2*mt) + 2 * np.outer(self.mrp,self.mrp)
        return 0.25 * np.transpose( B @ omega )
    
    # Method for obtaining the angle of rotation in degrees described by
    # this attitude coordinate.
    def get_angle(self):
        return 4 * np.rad2deg( np.arctan( (np.linalg.norm( self.mrp )) ) )
    
    # Method for obtaining the axis of rotation that this attitude describes.
    # Note that if the attitude describes a transformation from the B frame 
    # to the N frame, or vice versa, the axis of rotation is the same in both
    # frames and which is why it is an eigenvector of the rotation matrix.
    def get_axis(self):
        axis = np.array([ self.mrp[0], self.mrp[1], self.mrp[2] ])
        axis = axis / np.linalg.norm(axis)
        return axis
