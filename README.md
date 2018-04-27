# non_parametric
This is the R code for several common non-parametric methods (kernel est., mean regression, quantile regression, boostraps) with both practical applications on data and simulations

{::nomarkdown}

<!-- HTML CODE-->
 <html><head>
	<TITLE>RGL model</TITLE>
    </head>
    <body onload="rgl.start();"> 
    
    <div align="center">
<script>/*
 * Copyright (C) 2009 Apple Inc. All Rights Reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY APPLE INC. ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL APPLE INC. OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * Copyright (2016) Duncan Murdoch - fixed CanvasMatrix4.ortho,
 * cleaned up.
 */
/*
    CanvasMatrix4 class
    This class implements a 4x4 matrix. It has functions which
    duplicate the functionality of the OpenGL matrix stack and
    glut functions.
    IDL:
    [
        Constructor(in CanvasMatrix4 matrix),           // copy passed matrix into new CanvasMatrix4
        Constructor(in sequence<float> array)           // create new CanvasMatrix4 with 16 floats (row major)
        Constructor()                                   // create new CanvasMatrix4 with identity matrix
    ]
    interface CanvasMatrix4 {
        attribute float m11;
        attribute float m12;
        attribute float m13;
        attribute float m14;
        attribute float m21;
        attribute float m22;
        attribute float m23;
        attribute float m24;
        attribute float m31;
        attribute float m32;
        attribute float m33;
        attribute float m34;
        attribute float m41;
        attribute float m42;
        attribute float m43;
        attribute float m44;
        void load(in CanvasMatrix4 matrix);                 // copy the values from the passed matrix
        void load(in sequence<float> array);                // copy 16 floats into the matrix
        sequence<float> getAsArray();                       // return the matrix as an array of 16 floats
        WebGLFloatArray getAsCanvasFloatArray();           // return the matrix as a WebGLFloatArray with 16 values
        void makeIdentity();                                // replace the matrix with identity
        void transpose();                                   // replace the matrix with its transpose
        void invert();                                      // replace the matrix with its inverse
        void translate(in float x, in float y, in float z); // multiply the matrix by passed translation values on the right
        void scale(in float x, in float y, in float z);     // multiply the matrix by passed scale values on the right
        void rotate(in float angle,                         // multiply the matrix by passed rotation values on the right
                    in float x, in float y, in float z);    // (angle is in degrees)
        void multRight(in CanvasMatrix matrix);             // multiply the matrix by the passed matrix on the right
        void multLeft(in CanvasMatrix matrix);              // multiply the matrix by the passed matrix on the left
        void ortho(in float left, in float right,           // multiply the matrix by the passed ortho values on the right
                   in float bottom, in float top,
                   in float near, in float far);
        void frustum(in float left, in float right,         // multiply the matrix by the passed frustum values on the right
                     in float bottom, in float top,
                     in float near, in float far);
        void perspective(in float fovy, in float aspect,    // multiply the matrix by the passed perspective values on the right
                         in float zNear, in float zFar);
        void lookat(in float eyex, in float eyey, in float eyez,    // multiply the matrix by the passed lookat
                    in float ctrx, in float ctry, in float ctrz,    // values on the right
                    in float upx, in float upy, in float upz);
    }
*/
CanvasMatrix4 = function(m)
{
    if (typeof m == 'object') {
        if ("length" in m && m.length >= 16) {
            this.load(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
            return;
        }
        else if (m instanceof CanvasMatrix4) {
            this.load(m);
            return;
        }
    }
    this.makeIdentity();
};
CanvasMatrix4.prototype.load = function()
{
    if (arguments.length == 1 && typeof arguments[0] == 'object') {
        var matrix = arguments[0];
        if ("length" in matrix && matrix.length == 16) {
            this.m11 = matrix[0];
            this.m12 = matrix[1];
            this.m13 = matrix[2];
            this.m14 = matrix[3];
            this.m21 = matrix[4];
            this.m22 = matrix[5];
            this.m23 = matrix[6];
            this.m24 = matrix[7];
            this.m31 = matrix[8];
            this.m32 = matrix[9];
            this.m33 = matrix[10];
            this.m34 = matrix[11];
            this.m41 = matrix[12];
            this.m42 = matrix[13];
            this.m43 = matrix[14];
            this.m44 = matrix[15];
            return;
        }
        if (arguments[0] instanceof CanvasMatrix4) {
            this.m11 = matrix.m11;
            this.m12 = matrix.m12;
            this.m13 = matrix.m13;
            this.m14 = matrix.m14;
            this.m21 = matrix.m21;
            this.m22 = matrix.m22;
            this.m23 = matrix.m23;
            this.m24 = matrix.m24;
            this.m31 = matrix.m31;
            this.m32 = matrix.m32;
            this.m33 = matrix.m33;
            this.m34 = matrix.m34;
            this.m41 = matrix.m41;
            this.m42 = matrix.m42;
            this.m43 = matrix.m43;
            this.m44 = matrix.m44;
            return;
        }
    }
    this.makeIdentity();
};
CanvasMatrix4.prototype.getAsArray = function()
{
    return [
        this.m11, this.m12, this.m13, this.m14,
        this.m21, this.m22, this.m23, this.m24,
        this.m31, this.m32, this.m33, this.m34,
        this.m41, this.m42, this.m43, this.m44
    ];
};
CanvasMatrix4.prototype.getAsWebGLFloatArray = function()
{
    return new WebGLFloatArray(this.getAsArray());
};
CanvasMatrix4.prototype.makeIdentity = function()
{
    this.m11 = 1;
    this.m12 = 0;
    this.m13 = 0;
    this.m14 = 0;
    this.m21 = 0;
    this.m22 = 1;
    this.m23 = 0;
    this.m24 = 0;
    this.m31 = 0;
    this.m32 = 0;
    this.m33 = 1;
    this.m34 = 0;
    this.m41 = 0;
    this.m42 = 0;
    this.m43 = 0;
    this.m44 = 1;
};
CanvasMatrix4.prototype.transpose = function()
{
    var tmp = this.m12;
    this.m12 = this.m21;
    this.m21 = tmp;
    tmp = this.m13;
    this.m13 = this.m31;
    this.m31 = tmp;
    tmp = this.m14;
    this.m14 = this.m41;
    this.m41 = tmp;
    tmp = this.m23;
    this.m23 = this.m32;
    this.m32 = tmp;
    tmp = this.m24;
    this.m24 = this.m42;
    this.m42 = tmp;
    tmp = this.m34;
    this.m34 = this.m43;
    this.m43 = tmp;
};
CanvasMatrix4.prototype.invert = function()
{
    // Calculate the 4x4 determinant
    // If the determinant is zero,
    // then the inverse matrix is not unique.
    var det = this._determinant4x4();
    if (Math.abs(det) < 1e-8)
        return null;
    this._makeAdjoint();
    // Scale the adjoint matrix to get the inverse
    this.m11 /= det;
    this.m12 /= det;
    this.m13 /= det;
    this.m14 /= det;
    this.m21 /= det;
    this.m22 /= det;
    this.m23 /= det;
    this.m24 /= det;
    this.m31 /= det;
    this.m32 /= det;
    this.m33 /= det;
    this.m34 /= det;
    this.m41 /= det;
    this.m42 /= det;
    this.m43 /= det;
    this.m44 /= det;
};
CanvasMatrix4.prototype.translate = function(x,y,z)
{
    if (x === undefined)
        x = 0;
    if (y === undefined)
        y = 0;
    if (z === undefined)
        z = 0;
    var matrix = new CanvasMatrix4();
    matrix.m41 = x;
    matrix.m42 = y;
    matrix.m43 = z;
    this.multRight(matrix);
};
CanvasMatrix4.prototype.scale = function(x,y,z)
{
    if (x === undefined)
        x = 1;
    if (z === undefined) {
        if (y === undefined) {
            y = x;
            z = x;
        }
        else
            z = 1;
    }
    else if (y === undefined)
        y = x;
    var matrix = new CanvasMatrix4();
    matrix.m11 = x;
    matrix.m22 = y;
    matrix.m33 = z;
    this.multRight(matrix);
};
CanvasMatrix4.prototype.rotate = function(angle,x,y,z)
{
    // angles are in degrees. Switch to radians
    angle = angle / 180 * Math.PI;
    angle /= 2;
    var sinA = Math.sin(angle);
    var cosA = Math.cos(angle);
    var sinA2 = sinA * sinA;
    // normalize
    var length = Math.sqrt(x * x + y * y + z * z);
    if (length === 0) {
        // bad vector, just use something reasonable
        x = 0;
        y = 0;
        z = 1;
    } else if (length != 1) {
        x /= length;
        y /= length;
        z /= length;
    }
    var mat = new CanvasMatrix4();
    // optimize case where axis is along major axis
    if (x == 1 && y === 0 && z === 0) {
        mat.m11 = 1;
        mat.m12 = 0;
        mat.m13 = 0;
        mat.m21 = 0;
        mat.m22 = 1 - 2 * sinA2;
        mat.m23 = 2 * sinA * cosA;
        mat.m31 = 0;
        mat.m32 = -2 * sinA * cosA;
        mat.m33 = 1 - 2 * sinA2;
        mat.m14 = mat.m24 = mat.m34 = 0;
        mat.m41 = mat.m42 = mat.m43 = 0;
        mat.m44 = 1;
    } else if (x === 0 && y == 1 && z === 0) {
        mat.m11 = 1 - 2 * sinA2;
        mat.m12 = 0;
        mat.m13 = -2 * sinA * cosA;
        mat.m21 = 0;
        mat.m22 = 1;
        mat.m23 = 0;
        mat.m31 = 2 * sinA * cosA;
        mat.m32 = 0;
        mat.m33 = 1 - 2 * sinA2;
        mat.m14 = mat.m24 = mat.m34 = 0;
        mat.m41 = mat.m42 = mat.m43 = 0;
        mat.m44 = 1;
    } else if (x === 0 && y === 0 && z == 1) {
        mat.m11 = 1 - 2 * sinA2;
        mat.m12 = 2 * sinA * cosA;
        mat.m13 = 0;
        mat.m21 = -2 * sinA * cosA;
        mat.m22 = 1 - 2 * sinA2;
        mat.m23 = 0;
        mat.m31 = 0;
        mat.m32 = 0;
        mat.m33 = 1;
        mat.m14 = mat.m24 = mat.m34 = 0;
        mat.m41 = mat.m42 = mat.m43 = 0;
        mat.m44 = 1;
    } else {
        var x2 = x*x;
        var y2 = y*y;
        var z2 = z*z;
        mat.m11 = 1 - 2 * (y2 + z2) * sinA2;
        mat.m12 = 2 * (x * y * sinA2 + z * sinA * cosA);
        mat.m13 = 2 * (x * z * sinA2 - y * sinA * cosA);
        mat.m21 = 2 * (y * x * sinA2 - z * sinA * cosA);
        mat.m22 = 1 - 2 * (z2 + x2) * sinA2;
        mat.m23 = 2 * (y * z * sinA2 + x * sinA * cosA);
        mat.m31 = 2 * (z * x * sinA2 + y * sinA * cosA);
        mat.m32 = 2 * (z * y * sinA2 - x * sinA * cosA);
        mat.m33 = 1 - 2 * (x2 + y2) * sinA2;
        mat.m14 = mat.m24 = mat.m34 = 0;
        mat.m41 = mat.m42 = mat.m43 = 0;
        mat.m44 = 1;
    }
    this.multRight(mat);
};
CanvasMatrix4.prototype.multRight = function(mat)
{
    var m11 = (this.m11 * mat.m11 + this.m12 * mat.m21 +
               this.m13 * mat.m31 + this.m14 * mat.m41);
    var m12 = (this.m11 * mat.m12 + this.m12 * mat.m22 +
               this.m13 * mat.m32 + this.m14 * mat.m42);
    var m13 = (this.m11 * mat.m13 + this.m12 * mat.m23 +
               this.m13 * mat.m33 + this.m14 * mat.m43);
    var m14 = (this.m11 * mat.m14 + this.m12 * mat.m24 +
               this.m13 * mat.m34 + this.m14 * mat.m44);
    var m21 = (this.m21 * mat.m11 + this.m22 * mat.m21 +
               this.m23 * mat.m31 + this.m24 * mat.m41);
    var m22 = (this.m21 * mat.m12 + this.m22 * mat.m22 +
               this.m23 * mat.m32 + this.m24 * mat.m42);
    var m23 = (this.m21 * mat.m13 + this.m22 * mat.m23 +
               this.m23 * mat.m33 + this.m24 * mat.m43);
    var m24 = (this.m21 * mat.m14 + this.m22 * mat.m24 +
               this.m23 * mat.m34 + this.m24 * mat.m44);
    var m31 = (this.m31 * mat.m11 + this.m32 * mat.m21 +
               this.m33 * mat.m31 + this.m34 * mat.m41);
    var m32 = (this.m31 * mat.m12 + this.m32 * mat.m22 +
               this.m33 * mat.m32 + this.m34 * mat.m42);
    var m33 = (this.m31 * mat.m13 + this.m32 * mat.m23 +
               this.m33 * mat.m33 + this.m34 * mat.m43);
    var m34 = (this.m31 * mat.m14 + this.m32 * mat.m24 +
               this.m33 * mat.m34 + this.m34 * mat.m44);
    var m41 = (this.m41 * mat.m11 + this.m42 * mat.m21 +
               this.m43 * mat.m31 + this.m44 * mat.m41);
    var m42 = (this.m41 * mat.m12 + this.m42 * mat.m22 +
               this.m43 * mat.m32 + this.m44 * mat.m42);
    var m43 = (this.m41 * mat.m13 + this.m42 * mat.m23 +
               this.m43 * mat.m33 + this.m44 * mat.m43);
    var m44 = (this.m41 * mat.m14 + this.m42 * mat.m24 +
               this.m43 * mat.m34 + this.m44 * mat.m44);
    this.m11 = m11;
    this.m12 = m12;
    this.m13 = m13;
    this.m14 = m14;
    this.m21 = m21;
    this.m22 = m22;
    this.m23 = m23;
    this.m24 = m24;
    this.m31 = m31;
    this.m32 = m32;
    this.m33 = m33;
    this.m34 = m34;
    this.m41 = m41;
    this.m42 = m42;
    this.m43 = m43;
    this.m44 = m44;
};
CanvasMatrix4.prototype.multLeft = function(mat)
{
    var m11 = (mat.m11 * this.m11 + mat.m12 * this.m21 +
               mat.m13 * this.m31 + mat.m14 * this.m41);
    var m12 = (mat.m11 * this.m12 + mat.m12 * this.m22 +
               mat.m13 * this.m32 + mat.m14 * this.m42);
    var m13 = (mat.m11 * this.m13 + mat.m12 * this.m23 +
               mat.m13 * this.m33 + mat.m14 * this.m43);
    var m14 = (mat.m11 * this.m14 + mat.m12 * this.m24 +
               mat.m13 * this.m34 + mat.m14 * this.m44);
    var m21 = (mat.m21 * this.m11 + mat.m22 * this.m21 +
               mat.m23 * this.m31 + mat.m24 * this.m41);
    var m22 = (mat.m21 * this.m12 + mat.m22 * this.m22 +
               mat.m23 * this.m32 + mat.m24 * this.m42);
    var m23 = (mat.m21 * this.m13 + mat.m22 * this.m23 +
               mat.m23 * this.m33 + mat.m24 * this.m43);
    var m24 = (mat.m21 * this.m14 + mat.m22 * this.m24 +
               mat.m23 * this.m34 + mat.m24 * this.m44);
    var m31 = (mat.m31 * this.m11 + mat.m32 * this.m21 +
               mat.m33 * this.m31 + mat.m34 * this.m41);
    var m32 = (mat.m31 * this.m12 + mat.m32 * this.m22 +
               mat.m33 * this.m32 + mat.m34 * this.m42);
    var m33 = (mat.m31 * this.m13 + mat.m32 * this.m23 +
               mat.m33 * this.m33 + mat.m34 * this.m43);
    var m34 = (mat.m31 * this.m14 + mat.m32 * this.m24 +
               mat.m33 * this.m34 + mat.m34 * this.m44);
    var m41 = (mat.m41 * this.m11 + mat.m42 * this.m21 +
               mat.m43 * this.m31 + mat.m44 * this.m41);
    var m42 = (mat.m41 * this.m12 + mat.m42 * this.m22 +
               mat.m43 * this.m32 + mat.m44 * this.m42);
    var m43 = (mat.m41 * this.m13 + mat.m42 * this.m23 +
               mat.m43 * this.m33 + mat.m44 * this.m43);
    var m44 = (mat.m41 * this.m14 + mat.m42 * this.m24 +
               mat.m43 * this.m34 + mat.m44 * this.m44);
    this.m11 = m11;
    this.m12 = m12;
    this.m13 = m13;
    this.m14 = m14;
    this.m21 = m21;
    this.m22 = m22;
    this.m23 = m23;
    this.m24 = m24;
    this.m31 = m31;
    this.m32 = m32;
    this.m33 = m33;
    this.m34 = m34;
    this.m41 = m41;
    this.m42 = m42;
    this.m43 = m43;
    this.m44 = m44;
};
CanvasMatrix4.prototype.ortho = function(left, right, bottom, top, near, far)
{
    var tx = (left + right) / (left - right);
    var ty = (top + bottom) / (bottom - top);
    var tz = (far + near) / (near - far);
    var matrix = new CanvasMatrix4();
    matrix.m11 = 2 / (right - left);
    matrix.m12 = 0;
    matrix.m13 = 0;
    matrix.m14 = 0;
    matrix.m21 = 0;
    matrix.m22 = 2 / (top - bottom);
    matrix.m23 = 0;
    matrix.m24 = 0;
    matrix.m31 = 0;
    matrix.m32 = 0;
    matrix.m33 = -2 / (far - near);
    matrix.m34 = 0;
    matrix.m41 = tx;
    matrix.m42 = ty;
    matrix.m43 = tz;
    matrix.m44 = 1;
    this.multRight(matrix);
};
CanvasMatrix4.prototype.frustum = function(left, right, bottom, top, near, far)
{
    var matrix = new CanvasMatrix4();
    var A = (right + left) / (right - left);
    var B = (top + bottom) / (top - bottom);
    var C = -(far + near) / (far - near);
    var D = -(2 * far * near) / (far - near);
    matrix.m11 = (2 * near) / (right - left);
    matrix.m12 = 0;
    matrix.m13 = 0;
    matrix.m14 = 0;
    matrix.m21 = 0;
    matrix.m22 = 2 * near / (top - bottom);
    matrix.m23 = 0;
    matrix.m24 = 0;
    matrix.m31 = A;
    matrix.m32 = B;
    matrix.m33 = C;
    matrix.m34 = -1;
    matrix.m41 = 0;
    matrix.m42 = 0;
    matrix.m43 = D;
    matrix.m44 = 0;
    this.multRight(matrix);
};
CanvasMatrix4.prototype.perspective = function(fovy, aspect, zNear, zFar)
{
    var top = Math.tan(fovy * Math.PI / 360) * zNear;
    var bottom = -top;
    var left = aspect * bottom;
    var right = aspect * top;
    this.frustum(left, right, bottom, top, zNear, zFar);
};
CanvasMatrix4.prototype.lookat = function(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz)
{
    var matrix = new CanvasMatrix4();
    // Make rotation matrix
    // Z vector
    var zx = eyex - centerx;
    var zy = eyey - centery;
    var zz = eyez - centerz;
    var mag = Math.sqrt(zx * zx + zy * zy + zz * zz);
    if (mag) {
        zx /= mag;
        zy /= mag;
        zz /= mag;
    }
    // Y vector
    var yx = upx;
    var yy = upy;
    var yz = upz;
    // X vector = Y cross Z
    xx =  yy * zz - yz * zy;
    xy = -yx * zz + yz * zx;
    xz =  yx * zy - yy * zx;
    // Recompute Y = Z cross X
    yx = zy * xz - zz * xy;
    yy = -zx * xz + zz * xx;
    yx = zx * xy - zy * xx;
    // cross product gives area of parallelogram, which is < 1.0 for
    // non-perpendicular unit-length vectors; so normalize x, y here
    mag = Math.sqrt(xx * xx + xy * xy + xz * xz);
    if (mag) {
        xx /= mag;
        xy /= mag;
        xz /= mag;
    }
    mag = Math.sqrt(yx * yx + yy * yy + yz * yz);
    if (mag) {
        yx /= mag;
        yy /= mag;
        yz /= mag;
    }
    matrix.m11 = xx;
    matrix.m12 = xy;
    matrix.m13 = xz;
    matrix.m14 = 0;
    matrix.m21 = yx;
    matrix.m22 = yy;
    matrix.m23 = yz;
    matrix.m24 = 0;
    matrix.m31 = zx;
    matrix.m32 = zy;
    matrix.m33 = zz;
    matrix.m34 = 0;
    matrix.m41 = 0;
    matrix.m42 = 0;
    matrix.m43 = 0;
    matrix.m44 = 1;
    matrix.translate(-eyex, -eyey, -eyez);
    this.multRight(matrix);
};
// Support functions
CanvasMatrix4.prototype._determinant2x2 = function(a, b, c, d)
{
    return a * d - b * c;
};
CanvasMatrix4.prototype._determinant3x3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3)
{
    return a1 * this._determinant2x2(b2, b3, c2, c3) -
         b1 * this._determinant2x2(a2, a3, c2, c3) +
         c1 * this._determinant2x2(a2, a3, b2, b3);
};
CanvasMatrix4.prototype._determinant4x4 = function()
{
    var a1 = this.m11;
    var b1 = this.m12;
    var c1 = this.m13;
    var d1 = this.m14;
    var a2 = this.m21;
    var b2 = this.m22;
    var c2 = this.m23;
    var d2 = this.m24;
    var a3 = this.m31;
    var b3 = this.m32;
    var c3 = this.m33;
    var d3 = this.m34;
    var a4 = this.m41;
    var b4 = this.m42;
    var c4 = this.m43;
    var d4 = this.m44;
    return a1 * this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4) -
         b1 * this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4) +
         c1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4) -
         d1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
};
CanvasMatrix4.prototype._makeAdjoint = function()
{
    var a1 = this.m11;
    var b1 = this.m12;
    var c1 = this.m13;
    var d1 = this.m14;
    var a2 = this.m21;
    var b2 = this.m22;
    var c2 = this.m23;
    var d2 = this.m24;
    var a3 = this.m31;
    var b3 = this.m32;
    var c3 = this.m33;
    var d3 = this.m34;
    var a4 = this.m41;
    var b4 = this.m42;
    var c4 = this.m43;
    var d4 = this.m44;
    // Row column labeling reversed since we transpose rows & columns
    this.m11  =   this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
    this.m21  = - this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
    this.m31  =   this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
    this.m41  = - this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
    this.m12  = - this._determinant3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
    this.m22  =   this._determinant3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
    this.m32  = - this._determinant3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
    this.m42  =   this._determinant3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);
    this.m13  =   this._determinant3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
    this.m23  = - this._determinant3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
    this.m33  =   this._determinant3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
    this.m43  = - this._determinant3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);
    this.m14  = - this._determinant3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
    this.m24  =   this._determinant3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
    this.m34  = - this._determinant3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
    this.m44  =   this._determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
};</script>
<script>
rglwidgetClass = function() {
    this.canvas = null;
    this.userMatrix = new CanvasMatrix4();
    this.types = [];
    this.prMatrix = new CanvasMatrix4();
    this.mvMatrix = new CanvasMatrix4();
    this.vp = null;
    this.prmvMatrix = null;
    this.origs = null;
    this.gl = null;
    this.scene = null;
};
(function() {
    this.multMV = function(M, v) {
        return [ M.m11 * v[0] + M.m12 * v[1] + M.m13 * v[2] + M.m14 * v[3],
                 M.m21 * v[0] + M.m22 * v[1] + M.m23 * v[2] + M.m24 * v[3],
                 M.m31 * v[0] + M.m32 * v[1] + M.m33 * v[2] + M.m34 * v[3],
                 M.m41 * v[0] + M.m42 * v[1] + M.m43 * v[2] + M.m44 * v[3]
               ];
    };
    this.vlen = function(v) {
      return Math.sqrt(this.dotprod(v, v));
    };
    this.dotprod = function(a, b) {
      return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    this.xprod = function(a, b) {
      return [a[1]*b[2] - a[2]*b[1],
          a[2]*b[0] - a[0]*b[2],
          a[0]*b[1] - a[1]*b[0]];
    };
    this.cbind = function(a, b) {
      if (b.length < a.length)
        b = this.repeatToLen(b, a.length);
      else if (a.length < b.length)
        a = this.repeatToLen(a, b.length);
      return a.map(function(currentValue, index, array) {
            return currentValue.concat(b[index]);
      });
    };
    this.swap = function(a, i, j) {
      var temp = a[i];
      a[i] = a[j];
      a[j] = temp;
    };
    this.flatten = function(a) {
      return [].concat.apply([], a);
    };
    /* set element of 1d or 2d array as if it was flattened.  Column major, zero based! */
    this.setElement = function(a, i, value) {
      if (Array.isArray(a[0])) {
        var dim = a.length,
            col = Math.floor(i/dim),
            row = i % dim;
        a[row][col] = value;
      } else {
        a[i] = value;
      }
    };
    this.transpose = function(a) {
      var newArray = [],
          n = a.length,
          m = a[0].length,
          i;
      for(i = 0; i < m; i++){
        newArray.push([]);
      }
      for(i = 0; i < n; i++){
        for(var j = 0; j < m; j++){
          newArray[j].push(a[i][j]);
        }
      }
      return newArray;
    };
    this.sumsq = function(x) {
      var result = 0, i;
      for (i=0; i < x.length; i++)
        result += x[i]*x[i];
      return result;
    };
    this.toCanvasMatrix4 = function(mat) {
      if (mat instanceof CanvasMatrix4)
        return mat;
      var result = new CanvasMatrix4();
      mat = this.flatten(this.transpose(mat));
      result.load(mat);
      return result;
    };
    this.stringToRgb = function(s) {
      s = s.replace("#", "");
      var bigint = parseInt(s, 16);
      return [((bigint >> 16) & 255)/255,
              ((bigint >> 8) & 255)/255,
               (bigint & 255)/255];
    };
    this.componentProduct = function(x, y) {
      if (typeof y === "undefined") {
        this.alertOnce("Bad arg to componentProduct");
      }
      var result = new Float32Array(3), i;
      for (i = 0; i<3; i++)
        result[i] = x[i]*y[i];
      return result;
    };
    this.getPowerOfTwo = function(value) {
      var pow = 1;
      while(pow<value) {
        pow *= 2;
      }
      return pow;
    };
    this.unique = function(arr) {
      arr = [].concat(arr);
      return arr.filter(function(value, index, self) {
        return self.indexOf(value) === index;
      });
    };
    this.repeatToLen = function(arr, len) {
      arr = [].concat(arr);
      while (arr.length < len/2)
        arr = arr.concat(arr);
      return arr.concat(arr.slice(0, len - arr.length));
    };
    this.alertOnce = function(msg) {
      if (typeof this.alerted !== "undefined")
        return;
      this.alerted = true;
      alert(msg);
    };
    this.f_is_lit = 1;
    this.f_is_smooth = 2;
    this.f_has_texture = 4;
    this.f_is_indexed = 8;
    this.f_depth_sort = 16;
    this.f_fixed_quads = 32;
    this.f_is_transparent = 64;
    this.f_is_lines = 128;
    this.f_sprites_3d = 256;
    this.f_sprite_3d = 512;
    this.f_is_subscene = 1024;
    this.f_is_clipplanes = 2048;
    this.f_fixed_size = 4096;
    this.f_is_points = 8192;
    this.f_is_twosided = 16384;
    this.whichList = function(id) {
      var obj = this.getObj(id),
          flags = obj.flags;
        if (obj.type === "light")
          return "lights";
        if (flags & this.f_is_subscene)
            return "subscenes";
        if (flags & this.f_is_clipplanes)
            return "clipplanes";
        if (flags & this.f_is_transparent)
            return "transparent";
        return "opaque";
    };
    this.getObj = function(id) {
      if (typeof id !== "number") {
        this.alertOnce("getObj id is "+typeof id);
      }
      return this.scene.objects[id];
    };
    this.getIdsByType = function(type, subscene) {
      var
        result = [], i, self = this;
      if (typeof subscene === "undefined") {
        Object.keys(this.scene.objects).forEach(
          function(key) {
            key = parseInt(key, 10);
            if (self.getObj(key).type === type)
              result.push(key);
          });
      } else {
        ids = this.getObj(subscene).objects;
        for (i=0; i < ids.length; i++) {
          if (this.getObj(ids[i]).type === type) {
            result.push(ids[i]);
          }
        }
      }
      return result;
    };
    this.getMaterial = function(id, property) {
      var obj = this.getObj(id),
          mat = obj.material[property];
      if (typeof mat === "undefined")
          mat = this.scene.material[property];
      return mat;
    };
    this.inSubscene = function(id, subscene) {
      return this.getObj(subscene).objects.indexOf(id) > -1;
    };
    this.addToSubscene = function(id, subscene) {
      var thelist,
          thesub = this.getObj(subscene),
          ids = [id],
          obj = this.getObj(id), i;
      if (typeof obj.newIds !== "undefined") {
        ids = ids.concat(obj.newIds);
      }
      for (i = 0; i < ids.length; i++) {
        id = ids[i];
        if (thesub.objects.indexOf(id) == -1) {
          thelist = this.whichList(id);
          thesub.objects.push(id);
          thesub[thelist].push(id);
        }
      }
    };
    this.delFromSubscene = function(id, subscene) {
      var thelist,
          thesub = this.getObj(subscene),
          obj = this.getObj(id),
          ids = [id], i;
      if (typeof obj.newIds !== "undefined")
        ids = ids.concat(obj.newIds);
      for (j=0; j<ids.length;j++) {
        id = ids[j];
        i = thesub.objects.indexOf(id);
        if (i > -1) {
          thesub.objects.splice(i, 1);
          thelist = this.whichList(id);
          i = thesub[thelist].indexOf(id);
          thesub[thelist].splice(i, 1);
        }
      }
    };
    this.setSubsceneEntries = function(ids, subsceneid) {
      var sub = this.getObj(subsceneid);
      sub.objects = ids;
      this.initSubscene(subsceneid);
    };
    this.getSubsceneEntries = function(subscene) {
      return this.getObj(subscene).objects;
    };
    this.getChildSubscenes = function(subscene) {
      return this.getObj(subscene).subscenes;
    };
    this.startDrawing = function() {
    	var value = this.drawing;
    	this.drawing = true;
    	return value;
    };
    
    this.stopDrawing = function(saved) {
      this.drawing = saved;
      if (!saved && this.gl && this.gl.isContextLost())
        this.restartCanvas();
    };
    
    this.getVertexShader = function(id) {
      var obj = this.getObj(id),
          userShader = obj.userVertexShader,
          flags = obj.flags,
          type = obj.type,
          is_lit = flags & this.f_is_lit,
          has_texture = flags & this.f_has_texture,
          fixed_quads = flags & this.f_fixed_quads,
          sprites_3d = flags & this.f_sprites_3d,
          sprite_3d = flags & this.f_sprite_3d,
          nclipplanes = this.countClipplanes(),
          fixed_size = flags & this.f_fixed_size,
          is_points = flags & this.f_is_points,
          is_twosided = flags & this.f_is_twosided,
          result;
      if (type === "clipplanes" || sprites_3d) return;
      
      if (typeof userShader !== "undefined") return userShader;
      result = "  /* ****** "+type+" object "+id+" vertex shader ****** */\n"+
      "  attribute vec3 aPos;\n"+
      "  attribute vec4 aCol;\n"+
      " uniform mat4 mvMatrix;\n"+
      " uniform mat4 prMatrix;\n"+
      " varying vec4 vCol;\n"+
      " varying vec4 vPosition;\n";
      if ((is_lit && !fixed_quads) || sprite_3d)
        result = result + "  attribute vec3 aNorm;\n"+
                          " uniform mat4 normMatrix;\n"+
                          " varying vec3 vNormal;\n";
      if (has_texture || type === "text")
        result = result + " attribute vec2 aTexcoord;\n"+
                          " varying vec2 vTexcoord;\n";
      if (fixed_size)
        result = result + "  uniform vec2 textScale;\n";
      if (fixed_quads)
        result = result + "  attribute vec2 aOfs;\n";
      else if (sprite_3d)
        result = result + "  uniform vec3 uOrig;\n"+
                          "  uniform float uSize;\n"+
                          "  uniform mat4 usermat;\n";
                          
      if (is_twosided)
        result = result + "  attribute vec3 aPos1;\n"+
                          "  attribute vec3 aPos2;\n"+
                          "  varying float normz;\n";
      result = result + "  void main(void) {\n";
      if (nclipplanes || (!fixed_quads && !sprite_3d))
        result = result + "    vPosition = mvMatrix * vec4(aPos, 1.);\n";
      if (!fixed_quads && !sprite_3d)
        result = result + "    gl_Position = prMatrix * vPosition;\n";
      if (is_points) {
        var size = this.getMaterial(id, "size");
        result = result + "    gl_PointSize = "+size.toFixed(1)+";\n";
      }
      result = result + "    vCol = aCol;\n";
      if (is_lit && !fixed_quads && !sprite_3d)
        result = result + "    vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n";
      if (has_texture || type == "text")
        result = result + "    vTexcoord = aTexcoord;\n";
      if (fixed_size)
        result = result + "    vec4 pos = prMatrix * mvMatrix * vec4(aPos, 1.);\n"+
                          "   pos = pos/pos.w;\n"+
                          "   gl_Position = pos + vec4(aOfs*textScale, 0.,0.);\n";
      if (type == "sprites" && !fixed_size)
        result = result + "    vec4 pos = mvMatrix * vec4(aPos, 1.);\n"+
                          "   pos = pos/pos.w + vec4(aOfs, 0., 0.);\n"+
                          "   gl_Position = prMatrix*pos;\n";
      if (sprite_3d)
        result = result + "   vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n"+
                          "   vec4 pos = mvMatrix * vec4(uOrig, 1.);\n"+
                          "   vPosition = pos/pos.w + vec4(uSize*(vec4(aPos, 1.)*usermat).xyz,0.);\n"+
                          "   gl_Position = prMatrix * vPosition;\n";
      if (is_twosided)
        result = result + "   vec4 pos1 = prMatrix*(mvMatrix*vec4(aPos1, 1.));\n"+
                          "   pos1 = pos1/pos1.w - gl_Position/gl_Position.w;\n"+
                          "   vec4 pos2 = prMatrix*(mvMatrix*vec4(aPos2, 1.));\n"+
                          "   pos2 = pos2/pos2.w - gl_Position/gl_Position.w;\n"+
                          "   normz = pos1.x*pos2.y - pos1.y*pos2.x;\n";
      result = result + "  }\n";
      // console.log(result);
      return result;
    };
    this.getFragmentShader = function(id) {
      var obj = this.getObj(id),
          userShader = obj.userFragmentShader,
          flags = obj.flags,
          type = obj.type,
          is_lit = flags & this.f_is_lit,
          has_texture = flags & this.f_has_texture,
          fixed_quads = flags & this.f_fixed_quads,
          sprites_3d = flags & this.f_sprites_3d,
          is_twosided = (flags & this.f_is_twosided) > 0,
          nclipplanes = this.countClipplanes(), i,
          texture_format, nlights,
          result;
      if (type === "clipplanes" || sprites_3d) return;
      
      if (typeof userShader !== "undefined") return userShader;
      if (has_texture)
        texture_format = this.getMaterial(id, "textype");
      result = "/* ****** "+type+" object "+id+" fragment shader ****** */\n"+
               "#ifdef GL_ES\n"+
               "  precision highp float;\n"+
               "#endif\n"+
               "  varying vec4 vCol; // carries alpha\n"+
               "  varying vec4 vPosition;\n";
      if (has_texture || type === "text")
        result = result + "  varying vec2 vTexcoord;\n"+
                          " uniform sampler2D uSampler;\n";
      if (is_lit && !fixed_quads)
        result = result + "  varying vec3 vNormal;\n";
      for (i = 0; i < nclipplanes; i++)
        result = result + "  uniform vec4 vClipplane"+i+";\n";
      if (is_lit) {
        nlights = this.countLights();
        if (nlights)
            result = result + "  uniform mat4 mvMatrix;\n";
        else
            is_lit = false;
      }
      if (is_lit) {
        result = result + "   uniform vec3 emission;\n"+
                          "   uniform float shininess;\n";
        for (i=0; i < nlights; i++) {
          result = result + "   uniform vec3 ambient" + i + ";\n"+
                            "   uniform vec3 specular" + i +"; // light*material\n"+
                            "   uniform vec3 diffuse" + i + ";\n"+
                            "   uniform vec3 lightDir" + i + ";\n"+
                            "   uniform bool viewpoint" + i + ";\n"+
                            "   uniform bool finite" + i + ";\n";
        }
      }
      
      if (is_twosided) 
        result = result + "   uniform bool front;\n"+
                          "   varying float normz;\n";
      result = result + "  void main(void) {\n";
      for (i=0; i < nclipplanes;i++)
        result = result + "    if (dot(vPosition, vClipplane"+i+") < 0.0) discard;\n";
      if (fixed_quads) {
        result = result +   "    vec3 n = vec3(0., 0., 1.);\n";	
      } else if (is_lit) {
      	result = result +   "    vec3 n = normalize(vNormal);\n";
      }
      
      if (is_twosided) {
      	result = result +   "    if ((normz <= 0.) != front) discard;";
      }
      
      if (is_lit) {
        result = result + "    vec3 eye = normalize(-vPosition.xyz);\n"+
                          "   vec3 lightdir;\n"+
                          "   vec4 colDiff;\n"+
                          "   vec3 halfVec;\n"+
                          "   vec4 lighteffect = vec4(emission, 0.);\n"+
                          "   vec3 col;\n"+
                          "   float nDotL;\n";
        if (!fixed_quads) {
          result = result +   "   n = -faceforward(n, n, eye);\n";
        }
        for (i=0; i < nlights; i++) {
          result = result + "   colDiff = vec4(vCol.rgb * diffuse" + i + ", vCol.a);\n"+
                            "   lightdir = lightDir" + i + ";\n"+
                            "   if (!viewpoint" + i +")\n"+
                            "     lightdir = (mvMatrix * vec4(lightdir, 1.)).xyz;\n"+
                            "   if (!finite" + i + ") {\n"+
                            "     halfVec = normalize(lightdir + eye);\n"+
                            "   } else {\n"+
                            "     lightdir = normalize(lightdir - vPosition.xyz);\n"+
                            "     halfVec = normalize(lightdir + eye);\n"+
                            "   }\n"+
                            "    col = ambient" + i + ";\n"+
                            "   nDotL = dot(n, lightdir);\n"+
                            "   col = col + max(nDotL, 0.) * colDiff.rgb;\n"+
                            "   col = col + pow(max(dot(halfVec, n), 0.), shininess) * specular" + i + ";\n"+
                            "   lighteffect = lighteffect + vec4(col, colDiff.a);\n";
        }
      } else {
        result = result +   "   vec4 colDiff = vCol;\n"+
                            "    vec4 lighteffect = colDiff;\n";
      }
      if (type === "text")
        result = result +   "    vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n";
      if (has_texture) {
        result = result + {
            rgb:            "   vec4 textureColor = lighteffect*vec4(texture2D(uSampler, vTexcoord).rgb, 1.);\n",
            rgba:           "   vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n",
            alpha:          "   vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
                            "   float luminance = dot(vec3(1.,1.,1.), textureColor.rgb)/3.;\n"+
                            "   textureColor =  vec4(lighteffect.rgb, lighteffect.a*luminance);\n",
            luminance:      "   vec4 textureColor = vec4(lighteffect.rgb*dot(texture2D(uSampler, vTexcoord).rgb, vec3(1.,1.,1.))/3., lighteffect.a);\n",
          "luminance.alpha":"    vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
                            "   float luminance = dot(vec3(1.,1.,1.),textureColor.rgb)/3.;\n"+
                            "   textureColor = vec4(lighteffect.rgb*luminance, lighteffect.a*textureColor.a);\n"
          }[texture_format]+
                            "   gl_FragColor = textureColor;\n";
      } else if (type === "text") {
        result = result +   "    if (textureColor.a < 0.1)\n"+
                            "     discard;\n"+
                            "   else\n"+
                            "     gl_FragColor = textureColor;\n";
      } else
        result = result +   "   gl_FragColor = lighteffect;\n";
      result = result + "  }\n";
      // console.log(result);
      return result;
    };
    this.getShader = function(shaderType, code) {
        var gl = this.gl, shader;
        shader = gl.createShader(shaderType);
        gl.shaderSource(shader, code);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS) && !gl.isContextLost())
            alert(gl.getShaderInfoLog(shader));
        return shader;
    };
    this.handleLoadedTexture = function(texture, textureCanvas) { 
      var gl = this.gl || this.initGL();
      gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
      gl.bindTexture(gl.TEXTURE_2D, texture);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, textureCanvas);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
      gl.generateMipmap(gl.TEXTURE_2D);
      gl.bindTexture(gl.TEXTURE_2D, null);
    };
    this.loadImageToTexture = function(uri, texture) {
      var canvas = this.textureCanvas,
          ctx = canvas.getContext("2d"),
          image = new Image(),
          self = this;
       image.onload = function() {
         var w = image.width,
             h = image.height,
             canvasX = self.getPowerOfTwo(w),
             canvasY = self.getPowerOfTwo(h),
             gl = self.gl || self.initGL(),
             maxTexSize = gl.getParameter(gl.MAX_TEXTURE_SIZE);
         if (maxTexSize > 4096) maxTexSize = 4096;
         while (canvasX > 1 && canvasY > 1 && (canvasX > maxTexSize || canvasY > maxTexSize)) {
           canvasX /= 2;
           canvasY /= 2;
         }
         canvas.width = canvasX;
         canvas.height = canvasY;
         ctx.imageSmoothingEnabled = true;
         ctx.drawImage(image, 0, 0, canvasX, canvasY);
         self.handleLoadedTexture(texture, canvas);
         self.drawScene();
       };
       image.src = uri;
     };
    this.drawTextToCanvas = function(text, cex, family, font) {
       var canvasX, canvasY,
           textY,
           scaling = 20,
           textColour = "white",
           backgroundColour = "rgba(0,0,0,0)",
           canvas = this.textureCanvas,
           ctx = canvas.getContext("2d"),
           i, textHeights = [], widths = [], offset = 0, offsets = [],
           fontStrings = [],
           getFontString = function(i) {
             textHeights[i] = scaling*cex[i];
             var fontString = textHeights[i] + "px",
                 family0 = family[i],
                 font0 = font[i];
             if (family0 === "sans")
               family0 = "sans-serif";
             else if (family0 === "mono")
               family0 = "monospace";
             fontString = fontString + " " + family0;
             if (font0 === 2 || font0 === 4)
               fontString = "bold " + fontString;
             if (font0 === 3 || font0 === 4)
               fontString = "italic " + fontString;
             return fontString;
           };
       cex = this.repeatToLen(cex, text.length);
       family = this.repeatToLen(family, text.length);
       font = this.repeatToLen(font, text.length);
       canvasX = 1;
       for (i = 0; i < text.length; i++)  {
         ctx.font = fontStrings[i] = getFontString(i);
         widths[i] = ctx.measureText(text[i]).width;
         offset = offsets[i] = offset + 2*textHeights[i];
         canvasX = (widths[i] > canvasX) ? widths[i] : canvasX;
       }
       canvasX = this.getPowerOfTwo(canvasX);
       canvasY = this.getPowerOfTwo(offset);
       canvas.width = canvasX;
       canvas.height = canvasY;
       ctx.fillStyle = backgroundColour;
       ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);
       ctx.textBaseline = "alphabetic";
       for(i = 0; i < text.length; i++) {
         textY = offsets[i];
         ctx.font = fontStrings[i];
         ctx.fillStyle = textColour;
         ctx.textAlign = "left";
         ctx.fillText(text[i], 0,  textY);
       }
       return {canvasX:canvasX, canvasY:canvasY,
               widths:widths, textHeights:textHeights,
               offsets:offsets};
     };
    this.setViewport = function(id) {
       var gl = this.gl || this.initGL(),
         vp = this.getObj(id).par3d.viewport,
         x = vp.x*this.canvas.width,
         y = vp.y*this.canvas.height,
         width = vp.width*this.canvas.width,
         height = vp.height*this.canvas.height;
       this.vp = {x:x, y:y, width:width, height:height};
       gl.viewport(x, y, width, height);
       gl.scissor(x, y, width, height);
       gl.enable(gl.SCISSOR_TEST);
     };
    this.setprMatrix = function(id) {
       var subscene = this.getObj(id),
          embedding = subscene.embeddings.projection;
       if (embedding === "replace")
         this.prMatrix.makeIdentity();
       else
         this.setprMatrix(subscene.parent);
       if (embedding === "inherit")
         return;
       // This is based on the Frustum::enclose code from geom.cpp
       var bbox = subscene.par3d.bbox,
           scale = subscene.par3d.scale,
           ranges = [(bbox[1]-bbox[0])*scale[0]/2,
                     (bbox[3]-bbox[2])*scale[1]/2,
                     (bbox[5]-bbox[4])*scale[2]/2],
           radius = Math.sqrt(this.sumsq(ranges))*1.1; // A bit bigger to handle labels
       if (radius <= 0) radius = 1;
       var observer = subscene.par3d.observer,
           distance = observer[2],
           FOV = subscene.par3d.FOV, ortho = FOV === 0,
           t = ortho ? 1 : Math.tan(FOV*Math.PI/360),
           near = distance - radius,
           far = distance + radius,
           hlen,
           aspect = this.vp.width/this.vp.height,
           z = subscene.par3d.zoom;
       if (far < 0.)
         far = 1.;
       if (near < far/100.)
         near = far/100.;
       hlen = t*near;
       if (ortho) {
         if (aspect > 1)
           this.prMatrix.ortho(-hlen*aspect*z, hlen*aspect*z,
                          -hlen*z, hlen*z, near, far);
         else
           this.prMatrix.ortho(-hlen*z, hlen*z,
                          -hlen*z/aspect, hlen*z/aspect,
                          near, far);
       } else {
         if (aspect > 1)
           this.prMatrix.frustum(-hlen*aspect*z, hlen*aspect*z,
                          -hlen*z, hlen*z, near, far);
         else
           this.prMatrix.frustum(-hlen*z, hlen*z,
                          -hlen*z/aspect, hlen*z/aspect,
                          near, far);
       }
     };
    this.setmvMatrix = function(id) {
       var observer = this.getObj(id).par3d.observer;
       this.mvMatrix.makeIdentity();
       this.setmodelMatrix(id);
       this.mvMatrix.translate(-observer[0], -observer[1], -observer[2]);
     };
    this.setmodelMatrix = function(id) {
      var subscene = this.getObj(id),
          embedding = subscene.embeddings.model;
      if (embedding !== "inherit") {
        var scale = subscene.par3d.scale,
            bbox = subscene.par3d.bbox,
            center = [(bbox[0]+bbox[1])/2,
                      (bbox[2]+bbox[3])/2,
                      (bbox[4]+bbox[5])/2];
         this.mvMatrix.translate(-center[0], -center[1], -center[2]);
         this.mvMatrix.scale(scale[0], scale[1], scale[2]);
         this.mvMatrix.multRight( subscene.par3d.userMatrix );
       }
       if (embedding !== "replace")
         this.setmodelMatrix(subscene.parent);
     };
    this.setnormMatrix = function(subsceneid) {
       var self = this,
       recurse = function(id) {
         var sub = self.getObj(id),
             embedding = sub.embeddings.model;
         if (embedding !== "inherit") {
           var scale = sub.par3d.scale;
           self.normMatrix.scale(1/scale[0], 1/scale[1], 1/scale[2]);
           self.normMatrix.multRight(sub.par3d.userMatrix);
         }
         if (embedding !== "replace")
           recurse(sub.parent);
       };
       self.normMatrix.makeIdentity();
       recurse(subsceneid);
     };
    this.setprmvMatrix = function() {
       this.prmvMatrix = new CanvasMatrix4( this.mvMatrix );
       this.prmvMatrix.multRight( this.prMatrix );
     };
    this.countClipplanes = function() {
      return this.countObjs("clipplanes");
    };
    this.countLights = function() {
      return this.countObjs("light");
    };
    this.countObjs = function(type) {
      var self = this,
          bound = 0;
      Object.keys(this.scene.objects).forEach(
        function(key) {
          if (self.getObj(parseInt(key, 10)).type === type)
            bound = bound + 1;
        });
      return bound;
    };
    this.initSubscene = function(id) {
      var sub = this.getObj(id),
          i, obj;
      if (sub.type !== "subscene")
        return;
      sub.par3d.userMatrix = this.toCanvasMatrix4(sub.par3d.userMatrix);
      sub.par3d.listeners = [].concat(sub.par3d.listeners);
      sub.backgroundId = undefined;
      sub.subscenes = [];
      sub.clipplanes = [];
      sub.transparent = [];
      sub.opaque = [];
      sub.lights = [];
      for (i=0; i < sub.objects.length; i++) {
        obj = this.getObj(sub.objects[i]);
        if (typeof obj === "undefined") {
          sub.objects.splice(i, 1);
          i--;
        } else if (obj.type === "background")
          sub.backgroundId = obj.id;
        else
          sub[this.whichList(obj.id)].push(obj.id);
      }
    };
    this.copyObj = function(id, reuse) {
      var obj = this.getObj(id),
          prev = document.getElementById(reuse);
      if (prev !== null) {
        prev = prev.rglinstance;
        var
          prevobj = prev.getObj(id),
          fields = ["flags", "type",
                    "colors", "vertices", "centers",
                    "normals", "offsets",
                    "texts", "cex", "family", "font", "adj",
                    "material",
                    "radii",
                    "texcoords",
                    "userMatrix", "ids",
                    "dim",
                    "par3d", "userMatrix",
                    "viewpoint", "finite"],
          i;
        for (i = 0; i < fields.length; i++) {
          if (typeof prevobj[fields[i]] !== "undefined")
            obj[fields[i]] = prevobj[fields[i]];
        }
      } else
        console.warn("copyObj failed");
    };
    this.planeUpdateTriangles = function(id, bbox) {
      var perms = [[0,0,1], [1,2,2], [2,1,0]],
          x, xrow, elem, A, d, nhits, i, j, k, u, v, w, intersect, which, v0, v2, vx, reverse,
          face1 = [], face2 = [], normals = [],
          obj = this.getObj(id),
          nPlanes = obj.normals.length;
      obj.bbox = bbox;
      obj.vertices = [];
      obj.initialized = false;
      for (elem = 0; elem < nPlanes; elem++) {
//    Vertex Av = normal.getRecycled(elem);
        x = [];
        A = obj.normals[elem];
        d = obj.offsets[elem][0];
        nhits = 0;
        for (i=0; i<3; i++)
          for (j=0; j<2; j++)
            for (k=0; k<2; k++) {
              u = perms[0][i];
              v = perms[1][i];
              w = perms[2][i];
              if (A[w] !== 0.0) {
                intersect = -(d + A[u]*bbox[j+2*u] + A[v]*bbox[k+2*v])/A[w];
                if (bbox[2*w] < intersect && intersect < bbox[1+2*w]) {
                  xrow = [];
                  xrow[u] = bbox[j+2*u];
                  xrow[v] = bbox[k+2*v];
                  xrow[w] = intersect;
                  x.push(xrow);
                  face1[nhits] = j + 2*u;
                  face2[nhits] = k + 2*v;
                  nhits++;
                }
              }
            }
            if (nhits > 3) {
            /* Re-order the intersections so the triangles work */
              for (i=0; i<nhits-2; i++) {
                which = 0; /* initialize to suppress warning */
                for (j=i+1; j<nhits; j++) {
                  if (face1[i] == face1[j] || face1[i] == face2[j] ||
                      face2[i] == face1[j] || face2[i] == face2[j] ) {
                    which = j;
                    break;
                  }
                }
                if (which > i+1) {
                  this.swap(x, i+1, which);
                  this.swap(face1, i+1, which);
                  this.swap(face2, i+1, which);
                }
              }
            }
            if (nhits >= 3) {
      /* Put in order so that the normal points out the FRONT of the faces */
              v0 = [x[0][0] - x[1][0] , x[0][1] - x[1][1], x[0][2] - x[1][2]];
              v2 = [x[2][0] - x[1][0] , x[2][1] - x[1][1], x[2][2] - x[1][2]];
              /* cross-product */
              vx = this.xprod(v0, v2);
              reverse = this.dotprod(vx, A) > 0;
              for (i=0; i<nhits-2; i++) {
                obj.vertices.push(x[0]);
                normals.push(A);
                for (j=1; j<3; j++) {
                  obj.vertices.push(x[i + (reverse ? 3-j : j)]);
                  normals.push(A);
                }
              }
            }
      }
      obj.pnormals = normals;
    };
    this.initObj = function(id) {
      var obj = this.getObj(id),
          flags = obj.flags,
          type = obj.type,
          is_indexed = flags & this.f_is_indexed,
          is_lit = flags & this.f_is_lit,
          has_texture = flags & this.f_has_texture,
          fixed_quads = flags & this.f_fixed_quads,
          depth_sort = flags & this.f_depth_sort,
          sprites_3d = flags & this.f_sprites_3d,
          sprite_3d = flags & this.f_sprite_3d,
          fixed_size = flags & this.f_fixed_size,
          is_twosided = (flags & this.f_is_twosided) > 0,
          gl = this.gl || this.initGL(),
          texinfo, drawtype, nclipplanes, f, nrows,
          i,j,v,v1,v2, mat, uri, matobj, pass, pmode,
          dim, nx, nz, attr;
    if (typeof id !== "number") {
      this.alertOnce("initObj id is "+typeof id);
    }
    obj.initialized = true;
    if (type === "bboxdeco" || type === "subscene")
      return;
    if (type === "light") {
      obj.ambient = new Float32Array(obj.colors[0].slice(0,3));
      obj.diffuse = new Float32Array(obj.colors[1].slice(0,3));
      obj.specular = new Float32Array(obj.colors[2].slice(0,3));
      obj.lightDir = new Float32Array(obj.vertices[0]);
      return;
    }
    if (type === "clipplanes") {
      obj.vClipplane = this.flatten(this.cbind(obj.normals, obj.offsets));
      return;
    }
    if (type == "background" && typeof obj.ids !== "undefined") {
      obj.quad = this.flatten([].concat(obj.ids));
      return;
    }
    if (typeof obj.vertices === "undefined")
      obj.vertices = [];
    v = obj.vertices;
    obj.vertexCount = v.length;
    if (!obj.vertexCount) return;
    
    if (is_twosided) {
      if (typeof obj.userAttributes === "undefined")
        obj.userAttributes = {};
      v1 = Array(v.length);
      v2 = Array(v.length);  
      if (obj.type == "triangles" || obj.type == "quads") {
      	if (obj.type == "triangles")
      	  nrow = 3;
      	else
      	  nrow = 4;
        for (i=0; i<Math.floor(v.length/nrow); i++)
          for (j=0; j<nrow; j++) {
            v1[nrow*i + j] = v[nrow*i + ((j+1) % nrow)];
            v2[nrow*i + j] = v[nrow*i + ((j+2) % nrow)];
          }
      } else if (obj.type == "surface") {
        dim = obj.dim[0];
        nx = dim[0];
        nz = dim[1];
        for (j=0; j<nx; j++) {
          for (i=0; i<nz; i++) {
            if (i+1 < nz && j+1 < nx) { 
              v2[j + nx*i] = v[j + nx*(i+1)];
              v1[j + nx*i] = v[j+1 + nx*(i+1)];
            } else if (i+1 < nz) {
              v2[j + nx*i] = v[j-1 + nx*i];
              v1[j + nx*i] = v[j + nx*(i+1)];            	
            } else {
              v2[j + nx*i] = v[j + nx*(i-1)];
              v1[j + nx*i] = v[j-1 + nx*(i-1)];
            }
          }
        }
      }
      obj.userAttributes.aPos1 = v1;
      obj.userAttributes.aPos2 = v2;
    }
    if (!sprites_3d) {
      if (gl.isContextLost()) return;
      obj.prog = gl.createProgram();
      gl.attachShader(obj.prog, this.getShader( gl.VERTEX_SHADER,
        this.getVertexShader(id) ));
      gl.attachShader(obj.prog, this.getShader( gl.FRAGMENT_SHADER,
                      this.getFragmentShader(id) ));
      //  Force aPos to location 0, aCol to location 1
      gl.bindAttribLocation(obj.prog, 0, "aPos");
      gl.bindAttribLocation(obj.prog, 1, "aCol");
      gl.linkProgram(obj.prog);
      var linked = gl.getProgramParameter(obj.prog, gl.LINK_STATUS);
      if (!linked) {
        // An error occurred while linking
        var lastError = gl.getProgramInfoLog(obj.prog);
        console.warn("Error in program linking:" + lastError);
        gl.deleteProgram(obj.prog);
        return;
      }
    }
    if (type === "text") {
      texinfo = this.drawTextToCanvas(obj.texts,
                                      this.flatten(obj.cex),
                                      this.flatten(obj.family),
                                      this.flatten(obj.family));
    }
    if (fixed_quads && !sprites_3d) {
      obj.ofsLoc = gl.getAttribLocation(obj.prog, "aOfs");
    }
    if (sprite_3d) {
      obj.origLoc = gl.getUniformLocation(obj.prog, "uOrig");
      obj.sizeLoc = gl.getUniformLocation(obj.prog, "uSize");
      obj.usermatLoc = gl.getUniformLocation(obj.prog, "usermat");
    }
    if (has_texture || type == "text") {
      obj.texture = gl.createTexture();
      obj.texLoc = gl.getAttribLocation(obj.prog, "aTexcoord");
      obj.sampler = gl.getUniformLocation(obj.prog, "uSampler");
    }
    if (has_texture) {
      mat = obj.material;
      if (typeof mat.uri !== "undefined")
        uri = mat.uri;
      else if (typeof mat.uriElementId === "undefined") {
        matobj = this.getObj(mat.uriId);
        if (typeof matobj !== "undefined") {
          uri = matobj.material.uri;
        } else {
          uri = "";
        }
      } else
        uri = document.getElementById(mat.uriElementId).rglinstance.getObj(mat.uriId).material.uri;
      this.loadImageToTexture(uri, obj.texture);
    }
    if (type === "text") {
      this.handleLoadedTexture(obj.texture, this.textureCanvas);
    }
    
    var stride = 3, nc, cofs, nofs, radofs, oofs, tofs, vnew;
    nc = obj.colorCount = obj.colors.length;
    if (nc > 1) {
      cofs = stride;
      stride = stride + 4;
      v = this.cbind(v, obj.colors);
    } else {
      cofs = -1;
      obj.onecolor = this.flatten(obj.colors);
    }
    if (typeof obj.normals !== "undefined") {
      nofs = stride;
      stride = stride + 3;
      v = this.cbind(v, typeof obj.pnormals !== "undefined" ? obj.pnormals : obj.normals);
    } else
      nofs = -1;
    if (typeof obj.radii !== "undefined") {
      radofs = stride;
      stride = stride + 1;
      // FIXME:  always concat the radii?
      if (obj.radii.length === v.length) {
        v = this.cbind(v, obj.radii);
      } else if (obj.radii.length === 1) {
        v = v.map(function(row, i, arr) { return row.concat(obj.radii[0]);});
      }
    } else
      radofs = -1;
    if (type == "sprites" && !sprites_3d) {
      tofs = stride;
      stride += 2;
      oofs = stride;
      stride += 2;
      vnew = new Array(4*v.length);
      var rescale = fixed_size ? 72 : 1,
          size = obj.radii, s = rescale*size[0]/2;
      for (i=0; i < v.length; i++) {
        if (size.length > 1)
          s = rescale*size[i]/2;
        vnew[4*i]  = v[i].concat([0,0,-s,-s]);
        vnew[4*i+1]= v[i].concat([1,0, s,-s]);
        vnew[4*i+2]= v[i].concat([1,1, s, s]);
        vnew[4*i+3]= v[i].concat([0,1,-s, s]);
      }
      v = vnew;
      obj.vertexCount = v.length;
    } else if (type === "text") {
      tofs = stride;
      stride += 2;
      oofs = stride;
      stride += 2;
      vnew = new Array(4*v.length);
      for (i=0; i < v.length; i++) {
        vnew[4*i]  = v[i].concat([0,-0.5]).concat(obj.adj[0]);
        vnew[4*i+1]= v[i].concat([1,-0.5]).concat(obj.adj[0]);
        vnew[4*i+2]= v[i].concat([1, 1.5]).concat(obj.adj[0]);
        vnew[4*i+3]= v[i].concat([0, 1.5]).concat(obj.adj[0]);
        for (j=0; j < 4; j++) {
          v1 = vnew[4*i+j];
          v1[tofs+2] = 2*(v1[tofs]-v1[tofs+2])*texinfo.widths[i];
          v1[tofs+3] = 2*(v1[tofs+1]-v1[tofs+3])*texinfo.textHeights[i];
          v1[tofs] *= texinfo.widths[i]/texinfo.canvasX;
          v1[tofs+1] = 1.0-(texinfo.offsets[i] -
              v1[tofs+1]*texinfo.textHeights[i])/texinfo.canvasY;
          vnew[4*i+j] = v1;
        }
      }
      v = vnew;
      obj.vertexCount = v.length;
    } else if (typeof obj.texcoords !== "undefined") {
      tofs = stride;
      stride += 2;
      oofs = -1;
      v = this.cbind(v, obj.texcoords);
    } else {
      tofs = -1;
      oofs = -1;
    }
    
    if (typeof obj.userAttributes !== "undefined") {
      obj.userAttribOffsets = {};
      obj.userAttribLocations = {};
      obj.userAttribSizes = {};
      for (attr in obj.userAttributes) {
      	obj.userAttribLocations[attr] = gl.getAttribLocation(obj.prog, attr);
      	if (obj.userAttribLocations[attr] >= 0) { // Attribute may not have been used
      	  obj.userAttribOffsets[attr] = stride;
      	  v = this.cbind(v, obj.userAttributes[attr]);
      	  stride = v[0].length;
      	  obj.userAttribSizes[attr] = stride - obj.userAttribOffsets[attr];
      	}
      }
    }
    
    if (typeof obj.userUniforms !== "undefined") {
      obj.userUniformLocations = {};
      for (attr in obj.userUniforms) 
        obj.userUniformLocations[attr] = gl.getUniformLocation(obj.prog, attr);
    }
    if (stride !== v[0].length) {
      this.alertOnce("problem in stride calculation");
    }
    obj.vOffsets = {vofs:0, cofs:cofs, nofs:nofs, radofs:radofs, oofs:oofs, tofs:tofs, stride:stride};
    obj.values = new Float32Array(this.flatten(v));
    if (sprites_3d) {
      obj.userMatrix = new CanvasMatrix4(obj.userMatrix);
      obj.objects = this.flatten([].concat(obj.ids));
      is_lit = false;
    }
    if (is_lit && !fixed_quads) {
       obj.normLoc = gl.getAttribLocation(obj.prog, "aNorm");
    }
    nclipplanes = this.countClipplanes();
    if (nclipplanes && !sprites_3d) {
      obj.clipLoc = [];
      for (i=0; i < nclipplanes; i++)
        obj.clipLoc[i] = gl.getUniformLocation(obj.prog,"vClipplane" + i);
    }
    if (is_lit) {
      obj.emissionLoc = gl.getUniformLocation(obj.prog, "emission");
      obj.emission = new Float32Array(this.stringToRgb(this.getMaterial(id, "emission")));
      obj.shininessLoc = gl.getUniformLocation(obj.prog, "shininess");
      obj.shininess = this.getMaterial(id, "shininess");
      obj.nlights = this.countLights();
      obj.ambientLoc = [];
      obj.ambient = new Float32Array(this.stringToRgb(this.getMaterial(id, "ambient")));
      obj.specularLoc = [];
      obj.specular = new Float32Array(this.stringToRgb(this.getMaterial(id, "specular")));
      obj.diffuseLoc = [];
      obj.lightDirLoc = [];
      obj.viewpointLoc = [];
      obj.finiteLoc = [];
      for (i=0; i < obj.nlights; i++) {
        obj.ambientLoc[i] = gl.getUniformLocation(obj.prog, "ambient" + i);
        obj.specularLoc[i] = gl.getUniformLocation(obj.prog, "specular" + i);
        obj.diffuseLoc[i] = gl.getUniformLocation(obj.prog, "diffuse" + i);
        obj.lightDirLoc[i] = gl.getUniformLocation(obj.prog, "lightDir" + i);
        obj.viewpointLoc[i] = gl.getUniformLocation(obj.prog, "viewpoint" + i);
        obj.finiteLoc[i] = gl.getUniformLocation(obj.prog, "finite" + i);
      }
    }
    
    if (is_indexed) {
      obj.f = Array(2);
      for (pass = 0; pass < is_twosided + 1; pass++) {
      	if (type === "triangles" || type === "quads" || type === "surface")
      	  pmode = this.getMaterial(id, (pass === 0) ? "front" : "back");
      	else pmode = "filled";
      	if (pmode === "culled")
      	  continue;
        if (pmode === "points") {
      	  nrows = obj.vertexCount;
      	  f = Array(nrows);
      	  for (i=0; i < nrows; i++)
      	    f[i] = i;
        } else if ((type === "quads" || type === "text" ||
             type === "sprites") && !sprites_3d) {
          nrows = Math.floor(obj.vertexCount/4);
          if (pmode === "filled") {
            f = Array(6*nrows);
            for (i=0; i < nrows; i++) {
              f[6*i] = 4*i;
              f[6*i+1] = 4*i + 1;
              f[6*i+2] = 4*i + 2;
              f[6*i+3] = 4*i;
              f[6*i+4] = 4*i + 2;
              f[6*i+5] = 4*i + 3;
            }
          } else {
            f = Array(8*nrows);
            for (i=0; i < nrows; i++) {
              f[8*i] = 4*i;
              f[8*i+1] = 4*i + 1;
              f[8*i+2] = 4*i + 1;
              f[8*i+3] = 4*i + 2;
              f[8*i+4] = 4*i + 2;
              f[8*i+5] = 4*i + 3;
              f[8*i+6] = 4*i + 3;
              f[8*i+7] = 4*i;
            }
          }
        } else if (type === "triangles") {
          nrows = Math.floor(obj.vertexCount/3);
          if (pmode === "filled") {
            f = Array(3*nrows);
            for (i=0; i < f.length; i++) {
              f[i] = i;
            }
          } else if (pmode === "lines") {
            f = Array(6*nrows);
      	    for (i=0; i < nrows; i++) {
      	      f[6*i] = 3*i;
      	      f[6*i + 1] = 3*i + 1;
      	      f[6*i + 2] = 3*i + 1;
      	      f[6*i + 3] = 3*i + 2;
      	      f[6*i + 4] = 3*i + 2;  
      	      f[6*i + 5] = 3*i;      	  
      	    }
          }
        } else if (type === "spheres") {
          nrows = obj.vertexCount;
          f = Array(nrows);
          for (i=0; i < f.length; i++) {
            f[i] = i;
          }
        } else if (type === "surface") {
          dim = obj.dim[0];
          nx = dim[0];
          nz = dim[1];
          if (pmode === "filled") {
            f = [];
            for (j=0; j<nx-1; j++) {
              for (i=0; i<nz-1; i++) {
                f.push(j + nx*i,
                       j + nx*(i+1),
                       j + 1 + nx*(i+1),
                       j + nx*i,
                       j + 1 + nx*(i+1),
                       j + 1 + nx*i);
              } 
            }
          } else if (pmode === "lines") {
            f = [];
            for (j=0; j<nx; j++) {
              for (i=0; i<nz; i++) {
                if (i+1 < nz)
                  f.push(j + nx*i,
                         j + nx*(i+1));
                if (j+1 < nx)
                  f.push(j + nx*i,
                         j+1 + nx*i);
              }
            }
          }
        }
        obj.f[pass] = new Uint16Array(f);
        if (depth_sort) {
          drawtype = "DYNAMIC_DRAW";
        } else {
          drawtype = "STATIC_DRAW";
        }
      }
    }
    
    if (type !== "spheres" && !sprites_3d) {
      obj.buf = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
      gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW); //
    }
    if (is_indexed && type !== "spheres" && !sprites_3d) {
      obj.ibuf = Array(is_twosided + 1);
      obj.ibuf[0] = gl.createBuffer();
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[0]);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[0], gl[drawtype]);
      if (is_twosided) {
      	obj.ibuf[1] = gl.createBuffer();
      	gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[1]);
      	gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[1], gl[drawtype]);
      }
    }
    if (!sprites_3d) {
      obj.mvMatLoc = gl.getUniformLocation(obj.prog, "mvMatrix");
      obj.prMatLoc = gl.getUniformLocation(obj.prog, "prMatrix");
    }
    if (fixed_size) {
      obj.textScaleLoc = gl.getUniformLocation(obj.prog, "textScale");
    }
    if (is_lit && !sprites_3d) {
      obj.normMatLoc = gl.getUniformLocation(obj.prog, "normMatrix");
    }
    
    if (is_twosided) {
      obj.frontLoc = gl.getUniformLocation(obj.prog, "front");
    }
  };
    this.setDepthTest = function(id) {
      var gl = this.gl || this.initGL(),
          tests = {never: gl.NEVER,
                   less:  gl.LESS,
                   equal: gl.EQUAL,
                   lequal:gl.LEQUAL,
                   greater: gl.GREATER,
                   notequal: gl.NOTEQUAL,
                   gequal: gl.GEQUAL,
                   always: gl.ALWAYS},
           test = tests[this.getMaterial(id, "depth_test")];
      gl.depthFunc(test);
    };
    this.mode4type = {points : "POINTS",
                     linestrip : "LINE_STRIP",
                     abclines : "LINES",
                     lines : "LINES",
                     sprites : "TRIANGLES",
                     planes : "TRIANGLES",
                     text : "TRIANGLES",
                     quads : "TRIANGLES",
                     surface : "TRIANGLES",
                     triangles : "TRIANGLES"};
    this.drawObj = function(id, subsceneid) {
      var obj = this.getObj(id),
          subscene = this.getObj(subsceneid),
          flags = obj.flags,
          type = obj.type,
          is_indexed = flags & this.f_is_indexed,
          is_lit = flags & this.f_is_lit,
          has_texture = flags & this.f_has_texture,
          fixed_quads = flags & this.f_fixed_quads,
          depth_sort = flags & this.f_depth_sort,
          sprites_3d = flags & this.f_sprites_3d,
          sprite_3d = flags & this.f_sprite_3d,
          is_lines = flags & this.f_is_lines,
          is_points = flags & this.f_is_points,
          fixed_size = flags & this.f_fixed_size,
          is_twosided = (flags & this.f_is_twosided) > 0,
          gl = this.gl || this.initGL(),
          mat,
          sphereMV, baseofs, ofs, sscale, i, count, light,
          faces, pass, mode, pmode, attr,
          depthsort = function(i,j) { return depths[j] - depths[i]; };
      if (typeof id !== "number") {
        this.alertOnce("drawObj id is "+typeof id);
      }
      if (type === "planes") {
        if (obj.bbox !== subscene.par3d.bbox || !obj.initialized) {
          this.planeUpdateTriangles(id, subscene.par3d.bbox);
        }
      }
      if (!obj.initialized)
        this.initObj(id);
      if (type === "clipplanes") {
        count = obj.offsets.length;
        var IMVClip = [];
        for (i=0; i < count; i++) {
          IMVClip[i] = this.multMV(this.invMatrix, obj.vClipplane.slice(4*i, 4*(i+1)));
         }
         obj.IMVClip = IMVClip;
        return;
      }
      if (type === "light" || type === "bboxdeco" || !obj.vertexCount)
        return;
      this.setDepthTest(id);
      if (sprites_3d) {
        var norigs = obj.vertices.length,
            savenorm = new CanvasMatrix4(this.normMatrix);
        this.origs = obj.vertices;
        this.usermat = new Float32Array(obj.userMatrix.getAsArray());
        this.radii = obj.radii;
        this.normMatrix = subscene.spriteNormmat;
        for (this.iOrig=0; this.iOrig < norigs; this.iOrig++) {
          for (i=0; i < obj.objects.length; i++) {
            this.drawObj(obj.objects[i], subsceneid);
          }
        }
        this.normMatrix = savenorm;
        return;
      } else {
        gl.useProgram(obj.prog);
      }
      if (sprite_3d) {
        gl.uniform3fv(obj.origLoc, new Float32Array(this.origs[this.iOrig]));
        if (this.radii.length > 1) {
          gl.uniform1f(obj.sizeLoc, this.radii[this.iOrig][0]);
        } else {
          gl.uniform1f(obj.sizeLoc, this.radii[0][0]);
        }
        gl.uniformMatrix4fv(obj.usermatLoc, false, this.usermat);
      }
      
      if (type === "spheres") {
        gl.bindBuffer(gl.ARRAY_BUFFER, this.sphere.buf);
      } else {
        gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
      }
      gl.uniformMatrix4fv( obj.prMatLoc, false, new Float32Array(this.prMatrix.getAsArray()) );
      gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(this.mvMatrix.getAsArray()) );
      var clipcheck = 0,
          clipplaneids = subscene.clipplanes,
          clip, j;
      for (i=0; i < clipplaneids.length; i++) {
        clip = this.getObj(clipplaneids[i]);
        for (j=0; j < clip.offsets.length; j++) {
          gl.uniform4fv(obj.clipLoc[clipcheck + j], clip.IMVClip[j]);
        }
        clipcheck += clip.offsets.length;
      }
      if (typeof obj.clipLoc !== "undefined")
        for (i=clipcheck; i < obj.clipLoc.length; i++)
          gl.uniform4f(obj.clipLoc[i], 0,0,0,0);
      if (is_lit) {
        gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(this.normMatrix.getAsArray()) );
        gl.uniform3fv( obj.emissionLoc, obj.emission);
        gl.uniform1f( obj.shininessLoc, obj.shininess);
        for (i=0; i < subscene.lights.length; i++) {
          light = this.getObj(subscene.lights[i]);
          gl.uniform3fv( obj.ambientLoc[i], this.componentProduct(light.ambient, obj.ambient));
          gl.uniform3fv( obj.specularLoc[i], this.componentProduct(light.specular, obj.specular));
          gl.uniform3fv( obj.diffuseLoc[i], light.diffuse);
          gl.uniform3fv( obj.lightDirLoc[i], light.lightDir);
          gl.uniform1i( obj.viewpointLoc[i], light.viewpoint);
          gl.uniform1i( obj.finiteLoc[i], light.finite);
        }
        for (i=subscene.lights.length; i < obj.nlights; i++) {
          gl.uniform3f( obj.ambientLoc[i], 0,0,0);
          gl.uniform3f( obj.specularLoc[i], 0,0,0);
          gl.uniform3f( obj.diffuseLoc[i], 0,0,0);
        }
      }
      if (fixed_size) {
        gl.uniform2f( obj.textScaleLoc, 0.75/this.vp.width, 0.75/this.vp.height);
      }
      gl.enableVertexAttribArray( this.posLoc );
      var nc = obj.colorCount;
      count = obj.vertexCount;
      if (type === "spheres") {
        subscene = this.getObj(subsceneid);
        var scale = subscene.par3d.scale,
            scount = count;
        gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
        gl.enableVertexAttribArray(obj.normLoc );
        gl.vertexAttribPointer(obj.normLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
        gl.disableVertexAttribArray( this.colLoc );
        var sphereNorm = new CanvasMatrix4();
        sphereNorm.scale(scale[0], scale[1], scale[2]);
        sphereNorm.multRight(this.normMatrix);
        gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(sphereNorm.getAsArray()) );
        if (nc == 1) {
          gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
        }
        
        if (has_texture) {
          gl.enableVertexAttribArray( obj.texLoc );
          gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*this.sphere.vOffsets.stride, 
                                 4*this.sphere.vOffsets.tofs);
          gl.activeTexture(gl.TEXTURE0);
          gl.bindTexture(gl.TEXTURE_2D, obj.texture);
          gl.uniform1i( obj.sampler, 0);
        }
        for (i = 0; i < scount; i++) {
          sphereMV = new CanvasMatrix4();
          if (depth_sort) {
            baseofs = faces[i]*obj.vOffsets.stride;
          } else {
            baseofs = i*obj.vOffsets.stride;
          }
          ofs = baseofs + obj.vOffsets.radofs;
          sscale = obj.values[ofs];
          sphereMV.scale(sscale/scale[0], sscale/scale[1], sscale/scale[2]);
          sphereMV.translate(obj.values[baseofs],
                             obj.values[baseofs+1],
                             obj.values[baseofs+2]);
          sphereMV.multRight(this.mvMatrix);
          gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(sphereMV.getAsArray()) );
          if (nc > 1) {
            ofs = baseofs + obj.vOffsets.cofs;
            gl.vertexAttrib4f( this.colLoc, obj.values[ofs],
                                        obj.values[ofs+1],
                                       obj.values[ofs+2],
                                       obj.values[ofs+3] );
          }
          gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.sphere.ibuf);
          gl.drawElements(gl.TRIANGLES, this.sphere.sphereCount, gl.UNSIGNED_SHORT, 0);
        }
        return;
      } else {
        if (obj.colorCount === 1) {
          gl.disableVertexAttribArray( this.colLoc );
          gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
        } else {
          gl.enableVertexAttribArray( this.colLoc );
          gl.vertexAttribPointer(this.colLoc, 4, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.cofs);
        }
      }
      if (is_lit && obj.vOffsets.nofs > 0) {
        gl.enableVertexAttribArray( obj.normLoc );
        gl.vertexAttribPointer(obj.normLoc, 3, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.nofs);
      }
      if (has_texture || type === "text") {
        gl.enableVertexAttribArray( obj.texLoc );
        gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.tofs);
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, obj.texture);
        gl.uniform1i( obj.sampler, 0);
      }
      if (fixed_quads) {
        gl.enableVertexAttribArray( obj.ofsLoc );
        gl.vertexAttribPointer(obj.ofsLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.oofs);
      }
      
      if (typeof obj.userAttributes !== "undefined") {
      	for (attr in obj.userAttribSizes) {  // Not all attributes may have been used
      	  gl.enableVertexAttribArray( obj.userAttribLocations[attr] );
      	  gl.vertexAttribPointer( obj.userAttribLocations[attr], obj.userAttribSizes[attr],
      	  			  gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.userAttribOffsets[attr]);
      	}
      }
      
      if (typeof obj.userUniforms !== "undefined") {
      	for (attr in obj.userUniformLocations) {
      	  var loc = obj.userUniformLocations[attr];
      	  if (loc !== null) {
      	    var uniform = obj.userUniforms[attr];
      	    if (typeof uniform.length === "undefined")
      	      gl.uniform1f(loc, uniform);
      	    else if (typeof uniform[0].length === "undefined") {
      	      uniform = new Float32Array(uniform);
      	      switch(uniform.length) {
      	      	case 2: gl.uniform2fv(loc, uniform); break;
      	      	case 3: gl.uniform3fv(loc, uniform); break;
      	      	case 4: gl.uniform4fv(loc, uniform); break;
      	      	default: console.warn("bad uniform length");
      	      }
      	    } else if (uniform.length == 4 && uniform[0].length == 4) 
      	      gl.uniformMatrix4fv(loc, false, new Float32Array(uniform.getAsArray()));
      	    else
      	      console.warn("unsupported uniform matrix");
      	  }
      	}
      }
      for (pass = 0; pass < is_twosided + 1; pass++) {
      	if (type === "triangles" || type === "quads" || type === "surface")
      	  pmode = this.getMaterial(id, (pass === 0) ? "front" : "back");
      	else pmode = "filled";
        if (pmode === "culled")
          continue;
          
      	mode = this.mode4type[type];      
        if (depth_sort && pmode == "filled") {// Don't try depthsorting on wireframe or points
            var nfaces = obj.centers.length,
                z, w, frowsize;
            frowsize = Math.floor(obj.f[pass].length/nfaces);
            var depths = new Float32Array(nfaces);
            faces = new Array(nfaces);
            for(i=0; i<nfaces; i++) {
              z = this.prmvMatrix.m13*obj.centers[3*i] +
                  this.prmvMatrix.m23*obj.centers[3*i+1] +
                  this.prmvMatrix.m33*obj.centers[3*i+2] +
                  this.prmvMatrix.m43;
              w = this.prmvMatrix.m14*obj.centers[3*i] +
                  this.prmvMatrix.m24*obj.centers[3*i+1] +
                  this.prmvMatrix.m34*obj.centers[3*i+2] +
                  this.prmvMatrix.m44;
              depths[i] = z/w;
              faces[i] = i;
            }
            faces.sort(depthsort);
            if (type !== "spheres") {
              var f = new Uint16Array(obj.f[pass].length);
              for (i=0; i<nfaces; i++) {
                for (j=0; j<frowsize; j++) {
                  f[frowsize*i + j] = obj.f[pass][frowsize*faces[i] + j];
                }
              }
              gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
              gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, f, gl.DYNAMIC_DRAW);
            }
          }
      	
      	if (is_twosided)
      	  gl.uniform1i(obj.frontLoc, pass !== 0);
      	  
        if (is_indexed && type !== "spheres") {
          gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
        } else if (type === "spheres") {
          //  FIX ME!
        } 
      
        if (type === "sprites" || type === "text" || type === "quads") {
          count = count * 6/4;
        } else if (type === "surface") {
          count = obj.f[pass].length;
        }
      
        if (is_indexed) {
          count = obj.f[pass].length;
      	  if (pmode === "lines") {
      	    mode = "LINES";
      	    is_lines = true;
          } else if (pmode === "points") {
      	    mode = "POINTS";
          }
        }
        if (is_lines) {
          gl.lineWidth( this.getMaterial(id, "lwd") );
        }
        gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*obj.vOffsets.stride,  4*obj.vOffsets.vofs);
        if (is_indexed) {
          gl.drawElements(gl[mode], count, gl.UNSIGNED_SHORT, 0);
        } else {
          gl.drawArrays(gl[mode], 0, count);
        }
     }
   };
    this.drawBackground = function(id, subsceneid) {
      var gl = this.gl || this.initGL(),
          obj = this.getObj(id),
          bg, i;
      if (!obj.initialized)
        this.initObj(id);
      if (obj.colors.length) {
        bg = obj.colors[0];
        gl.clearColor(bg[0], bg[1], bg[2], bg[3]);
        gl.depthMask(true);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      }
      if (typeof obj.quad !== "undefined") {
        this.prMatrix.makeIdentity();
        this.mvMatrix.makeIdentity();
        gl.disable(gl.BLEND);
        gl.disable(gl.DEPTH_TEST);
        gl.depthMask(false);
        for (i=0; i < obj.quad.length; i++)
          this.drawObj(obj.quad[i], subsceneid);
      }
    };
    this.drawSubscene = function(subsceneid) {
      var gl = this.gl || this.initGL(),
          obj = this.getObj(subsceneid),
          objects = this.scene.objects,
          subids = obj.objects,
          subscene_has_faces = false,
          subscene_needs_sorting = false,
          flags, i;
      if (obj.par3d.skipRedraw)
        return;
      for (i=0; i < subids.length; i++) {
        flags = objects[subids[i]].flags;
        if (typeof flags !== "undefined") {
          subscene_has_faces |= (flags & this.f_is_lit)
                           & !(flags & this.f_fixed_quads);
          subscene_needs_sorting |= (flags & this.f_depth_sort);
        }
      }
      this.setViewport(subsceneid);
      if (typeof obj.backgroundId !== "undefined")
          this.drawBackground(obj.backgroundId, subsceneid);
      if (subids.length) {
        this.setprMatrix(subsceneid);
        this.setmvMatrix(subsceneid);
        if (subscene_has_faces) {
          this.setnormMatrix(subsceneid);
          if ((obj.flags & this.f_sprites_3d) &&
              typeof obj.spriteNormmat === "undefined") {
            obj.spriteNormmat = new CanvasMatrix4(this.normMatrix);
          }
        }
        if (subscene_needs_sorting)
          this.setprmvMatrix();
        gl.enable(gl.DEPTH_TEST);
        gl.depthMask(true);
        gl.disable(gl.BLEND);
        var clipids = obj.clipplanes;
        if (typeof clipids === "undefined") {
          console.warn("bad clipids");
        }
        if (clipids.length > 0) {
          this.invMatrix = new CanvasMatrix4(this.mvMatrix);
          this.invMatrix.invert();
          for (i = 0; i < clipids.length; i++)
            this.drawObj(clipids[i], subsceneid);
        }
        subids = obj.opaque;
        if (subids.length > 0) {
          for (i = 0; i < subids.length; i++) {
            this.drawObj(subids[i], subsceneid);
          }
        }
        subids = obj.transparent;
        if (subids.length > 0) {
          gl.depthMask(false);
          gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA,
                               gl.ONE, gl.ONE);
          gl.enable(gl.BLEND);
          for (i = 0; i < subids.length; i++) {
            this.drawObj(subids[i], subsceneid);
          }
        }
        subids = obj.subscenes;
        for (i = 0; i < subids.length; i++) {
          this.drawSubscene(subids[i]);
        }
      }
    };
    this.relMouseCoords = function(event) {
      var totalOffsetX = 0,
      totalOffsetY = 0,
      currentElement = this.canvas;
      do {
        totalOffsetX += currentElement.offsetLeft;
        totalOffsetY += currentElement.offsetTop;
        currentElement = currentElement.offsetParent;
      }
      while(currentElement);
      var canvasX = event.pageX - totalOffsetX,
          canvasY = event.pageY - totalOffsetY;
      return {x:canvasX, y:canvasY};
    };
    this.setMouseHandlers = function() {
      var self = this, activeSubscene, handler,
          handlers = {}, drag = 0;
      handlers.rotBase = 0;
      this.screenToVector = function(x, y) {
        var viewport = this.getObj(activeSubscene).par3d.viewport,
          width = viewport.width*this.canvas.width,
          height = viewport.height*this.canvas.height,
          radius = Math.max(width, height)/2.0,
          cx = width/2.0,
          cy = height/2.0,
          px = (x-cx)/radius,
          py = (y-cy)/radius,
          plen = Math.sqrt(px*px+py*py);
        if (plen > 1.e-6) {
          px = px/plen;
          py = py/plen;
        }
        var angle = (Math.SQRT2 - plen)/Math.SQRT2*Math.PI/2,
          z = Math.sin(angle),
          zlen = Math.sqrt(1.0 - z*z);
        px = px * zlen;
        py = py * zlen;
        return [px, py, z];
      };
      handlers.trackballdown = function(x,y) {
        var activeSub = this.getObj(activeSubscene),
            activeModel = this.getObj(this.useid(activeSub.id, "model")),
            i, l = activeModel.par3d.listeners;
        handlers.rotBase = this.screenToVector(x, y);
        this.saveMat = [];
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
        }
      };
      handlers.trackballmove = function(x,y) {
        var rotCurrent = this.screenToVector(x,y),
            rotBase = handlers.rotBase,
            dot = rotBase[0]*rotCurrent[0] +
                  rotBase[1]*rotCurrent[1] +
                  rotBase[2]*rotCurrent[2],
            angle = Math.acos( dot/this.vlen(rotBase)/this.vlen(rotCurrent) )*180.0/Math.PI,
            axis = this.xprod(rotBase, rotCurrent),
            objects = this.scene.objects,
            activeSub = this.getObj(activeSubscene),
            activeModel = this.getObj(this.useid(activeSub.id, "model")),
            l = activeModel.par3d.listeners,
            i;
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.par3d.userMatrix.load(objects[l[i]].saveMat);
          activeSub.par3d.userMatrix.rotate(angle, axis[0], axis[1], axis[2]);
        }
        this.drawScene();
      };
      handlers.trackballend = 0;
      handlers.axisdown = function(x,y) {
        handlers.rotBase = this.screenToVector(x, this.canvas.height/2);
        var activeSub = this.getObj(activeSubscene),
            activeModel = this.getObj(this.useid(activeSub.id, "model")),
            i, l = activeModel.par3d.listeners;
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
        }
      };
      handlers.axismove = function(x,y) {
        var rotCurrent = this.screenToVector(x, this.canvas.height/2),
            rotBase = handlers.rotBase,
            angle = (rotCurrent[0] - rotBase[0])*180/Math.PI,
            rotMat = new CanvasMatrix4();
        rotMat.rotate(angle, handlers.axis[0], handlers.axis[1], handlers.axis[2]);
        var activeSub = this.getObj(activeSubscene),
            activeModel = this.getObj(this.useid(activeSub.id, "model")),
            i, l = activeModel.par3d.listeners;
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.par3d.userMatrix.load(activeSub.saveMat);
          activeSub.par3d.userMatrix.multLeft(rotMat);
        }
        this.drawScene();
      };
      handlers.axisend = 0;
      handlers.y0zoom = 0;
      handlers.zoom0 = 0;
      handlers.zoomdown = function(x, y) {
        var activeSub = this.getObj(activeSubscene),
          activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
          i, l = activeProjection.par3d.listeners;
        handlers.y0zoom = y;
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.zoom0 = Math.log(activeSub.par3d.zoom);
        }
      };
      handlers.zoommove = function(x, y) {
        var activeSub = this.getObj(activeSubscene),
            activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
            i, l = activeProjection.par3d.listeners;
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.par3d.zoom = Math.exp(activeSub.zoom0 + (y-handlers.y0zoom)/this.canvas.height);
        }
        this.drawScene();
      };
      handlers.zoomend = 0;
      handlers.y0fov = 0;
      handlers.fovdown = function(x, y) {
        handlers.y0fov = y;
        var activeSub = this.getObj(activeSubscene),
          activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
          i, l = activeProjection.par3d.listeners;
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.fov0 = activeSub.par3d.FOV;
        }
      };
      handlers.fovmove = function(x, y) {
        var activeSub = this.getObj(activeSubscene),
            activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
            i, l = activeProjection.par3d.listeners;
        for (i = 0; i < l.length; i++) {
          activeSub = this.getObj(l[i]);
          activeSub.par3d.FOV = Math.max(1, Math.min(179, activeSub.fov0 +
             180*(y-handlers.y0fov)/this.canvas.height));
        }
        this.drawScene();
      };
      handlers.fovend = 0;
      this.canvas.onmousedown = function ( ev ){
        if (!ev.which) // Use w3c defns in preference to MS
        switch (ev.button) {
          case 0: ev.which = 1; break;
          case 1:
          case 4: ev.which = 2; break;
          case 2: ev.which = 3;
        }
        drag = ["left", "middle", "right"][ev.which-1];
        var coords = self.relMouseCoords(ev);
        coords.y = self.canvas.height-coords.y;
        activeSubscene = self.whichSubscene(coords);
        var sub = self.getObj(activeSubscene), f;
        handler = sub.par3d.mouseMode[drag];
        switch (handler) {
        case "xAxis":
          handler = "axis";
          handlers.axis = [1.0, 0.0, 0.0];
          break;
        case "yAxis":
          handler = "axis";
          handlers.axis = [0.0, 1.0, 0.0];
          break;
        case "zAxis":
          handler = "axis";
          handlers.axis = [0.0, 0.0, 1.0];
          break;
        }
        f = handlers[handler + "down"];
        if (f) {
          coords = self.translateCoords(activeSubscene, coords);
          f.call(self, coords.x, coords.y);
          ev.preventDefault();
        }
      };
      this.canvas.onmouseup = function ( ev ){
        if ( drag === 0 ) return;
        var f = handlers[handler + "up"];
        if (f)
          f();
        drag = 0;
      };
      this.canvas.onmouseout = this.canvas.onmouseup;
      this.canvas.onmousemove = function ( ev ) {
        if ( drag === 0 ) return;
        var f = handlers[handler + "move"];
        if (f) {
          var coords = self.relMouseCoords(ev);
          coords.y = self.canvas.height - coords.y;
          coords = self.translateCoords(activeSubscene, coords);
          f.call(self, coords.x, coords.y);
        }
      };
      handlers.wheelHandler = function(ev) {
        var del = 1.02, i;
        if (ev.shiftKey) del = 1.002;
        var ds = ((ev.detail || ev.wheelDelta) > 0) ? del : (1 / del);
        if (typeof activeSubscene === "undefined")
          activeSubscene = self.scene.rootSubscene;
        var activeSub = self.getObj(activeSubscene),
            activeProjection = self.getObj(self.useid(activeSub.id, "projection")),
            l = activeProjection.par3d.listeners;
        for (i = 0; i < l.length; i++) {
          activeSub = self.getObj(l[i]);
          activeSub.par3d.zoom *= ds;
        }
        self.drawScene();
        ev.preventDefault();
      };
      this.canvas.addEventListener("DOMMouseScroll", handlers.wheelHandler, false);
      this.canvas.addEventListener("mousewheel", handlers.wheelHandler, false);
    };
    this.useid = function(subsceneid, type) {
      var sub = this.getObj(subsceneid);
      if (sub.embeddings[type] === "inherit")
        return(this.useid(sub.parent, type));
      else
        return subsceneid;
    };
    this.inViewport = function(coords, subsceneid) {
      var viewport = this.getObj(subsceneid).par3d.viewport,
        x0 = coords.x - viewport.x*this.canvas.width,
        y0 = coords.y - viewport.y*this.canvas.height;
      return 0 <= x0 && x0 <= viewport.width*this.canvas.width &&
             0 <= y0 && y0 <= viewport.height*this.canvas.height;
    };
    this.whichSubscene = function(coords) {
      var self = this,
          recurse = function(subsceneid) {
            var subscenes = self.getChildSubscenes(subsceneid), i, id;
            for (i=0; i < subscenes.length; i++) {
              id = recurse(subscenes[i]);
              if (typeof(id) !== "undefined")
                return(id);
            }
            if (self.inViewport(coords, subsceneid))
              return(subsceneid);
            else
              return undefined;
          },
          rootid = this.scene.rootSubscene,
          result = recurse(rootid);
      if (typeof(result) === "undefined")
        result = rootid;
      return result;
    };
    this.translateCoords = function(subsceneid, coords) {
      var viewport = this.getObj(subsceneid).par3d.viewport;
      return {x: coords.x - viewport.x*this.canvas.width,
              y: coords.y - viewport.y*this.canvas.height};
    };
    this.initSphere = function() {
      var verts = this.scene.sphereVerts, 
          reuse = verts.reuse, result;
      if (typeof reuse !== "undefined") {
        var prev = document.getElementById(reuse).rglinstance.sphere;
        result = {values: prev.values, vOffsets: prev.vOffsets, it: prev.it};
      } else 
        result = {values: new Float32Array(this.flatten(this.cbind(this.transpose(verts.vb),
                    this.transpose(verts.texcoords)))),
                  it: new Uint16Array(this.flatten(this.transpose(verts.it))),
                  vOffsets: {vofs:0, cofs:-1, nofs:-1, radofs:-1, oofs:-1, 
                    tofs:3, stride:5}};
      result.sphereCount = result.it.length;
      this.sphere = result;
    };
    
    this.initSphereGL = function() {
      var gl = this.gl || this.initGL(), sphere = this.sphere;
      if (gl.isContextLost()) return;
      sphere.buf = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, sphere.buf);
      gl.bufferData(gl.ARRAY_BUFFER, sphere.values, gl.STATIC_DRAW);
      sphere.ibuf = gl.createBuffer();
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sphere.ibuf);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, sphere.it, gl.STATIC_DRAW);
      return;
    };
    this.initialize = function(el, x) {
      this.textureCanvas = document.createElement("canvas");
      this.textureCanvas.style.display = "block";
      this.scene = x;
      this.normMatrix = new CanvasMatrix4();
      this.saveMat = {};
      this.distance = null;
      this.posLoc = 0;
      this.colLoc = 1;
      if (el) {
        el.rglinstance = this;
        this.el = el;
        this.webGLoptions = el.rglinstance.scene.webGLoptions;
        this.initCanvas();
      }
    };
    this.restartCanvas = function() {
      var newcanvas = document.createElement("canvas");
      newcanvas.width = this.el.width;
      newcanvas.height = this.el.height;
      newcanvas.addEventListener("webglcontextrestored",
        this.onContextRestored, false);
      newcanvas.addEventListener("webglcontextlost",
        this.onContextLost, false);            
      while (this.el.firstChild) {
        this.el.removeChild(this.el.firstChild);
      }
      this.el.appendChild(newcanvas);
      this.canvas = newcanvas;
      this.gl = null;      	
    };
      
    this.initCanvas = function() {
      this.restartCanvas();
      var objs = this.scene.objects,
          self = this;
      Object.keys(objs).forEach(function(key){
        var id = parseInt(key, 10),
            obj = self.getObj(id);
        if (typeof obj.reuse !== "undefined")
          self.copyObj(id, obj.reuse);
      });
      Object.keys(objs).forEach(function(key){
        self.initSubscene(parseInt(key, 10));
      });
      this.setMouseHandlers();      
      this.initSphere();
      
      this.onContextRestored = function(event) {
        self.initGL();
        self.drawScene();
        // console.log("restored context for "+self.scene.rootSubscene);
      };
      
      this.onContextLost = function(event) {
        if (!self.drawing)
          self.restartCanvas();
        event.preventDefault();
      };
      
      this.initGL0();
      lazyLoadScene = function() {
      	if (self.isInBrowserViewport()) {
      	  if (!self.gl) {
      	    self.initGL();
      	  }
      	  self.drawScene();
      	}
      };
      window.addEventListener("DOMContentLoaded", lazyLoadScene, false);
      window.addEventListener("load", lazyLoadScene, false);
      window.addEventListener("resize", lazyLoadScene, false);
      window.addEventListener("scroll", lazyLoadScene, false);
    };
    /* this is only used by writeWebGL; rglwidget has
       no debug element and does the drawing in rglwidget.js */
    this.start = function() {
      if (typeof this.prefix !== "undefined") {
        this.debugelement = document.getElementById(this.prefix + "debug");
        this.debug("");
      }
      this.drag = 0;
      this.drawScene();
    };
    this.debug = function(msg, img) {
      if (typeof this.debugelement !== "undefined" && this.debugelement !== null) {
        this.debugelement.innerHTML = msg;
        if (typeof img !== "undefined") {
          this.debugelement.insertBefore(img, this.debugelement.firstChild);
        }
      } else if (msg !== "")
        alert(msg);
    };
    this.getSnapshot = function() {
      var img;
      if (typeof this.scene.snapshot !== "undefined") {
        img = document.createElement("img");
        img.src = this.scene.snapshot;
        img.alt = "Snapshot";
      }
      return img;
    };
    this.initGL0 = function() {
      if (!window.WebGLRenderingContext){
        alert("Your browser does not support WebGL. See http://get.webgl.org");
        return;
      }
    };
    
    this.isInBrowserViewport = function() {
      var rect = this.canvas.getBoundingClientRect(),
          windHeight = (window.innerHeight || document.documentElement.clientHeight),
          windWidth = (window.innerWidth || document.documentElement.clientWidth);
      return (
      	rect.top >= -windHeight &&
      	rect.left >= -windWidth &&
      	rect.bottom <= 2*windHeight &&
      	rect.right <= 2*windWidth);
    };
    this.initGL = function() {
      var self = this;
      if (this.gl) {
      	if (!this.drawing && this.gl.isContextLost())
          this.restartCanvas();
        else
          return this.gl;
      }
      // if (!this.isInBrowserViewport()) return; Return what??? At this point we know this.gl is null.
      this.canvas.addEventListener("webglcontextrestored",
        this.onContextRestored, false);
      this.canvas.addEventListener("webglcontextlost",
        this.onContextLost, false);      
      this.gl = this.canvas.getContext("webgl", this.webGLoptions) ||
               this.canvas.getContext("experimental-webgl", this.webGLoptions);
      var save = this.startDrawing();
      this.initSphereGL(); 
      Object.keys(this.scene.objects).forEach(function(key){
        self.initObj(parseInt(key, 10));
        });
      this.stopDrawing(save);
      return this.gl;
    };
    this.resize = function(el) {
      this.canvas.width = el.width;
      this.canvas.height = el.height;
    };
    this.drawScene = function() {
      var gl = this.gl || this.initGL(),
          save = this.startDrawing();
      gl.enable(gl.DEPTH_TEST);
      gl.depthFunc(gl.LEQUAL);
      gl.clearDepth(1.0);
      gl.clearColor(1,1,1,1);
      gl.depthMask(true); // Must be true before clearing depth buffer
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      this.drawSubscene(this.scene.rootSubscene);
      this.stopDrawing(save);
    };
    this.subsetSetter = function(el, control) {
      if (typeof control.subscenes === "undefined" ||
          control.subscenes === null)
        control.subscenes = this.scene.rootSubscene;
      var value = Math.round(control.value),
          subscenes = [].concat(control.subscenes),
          fullset = [].concat(control.fullset),
          i, j, entries, subsceneid, 
          adds = [], deletes = [],
          ismissing = function(x) {
            return fullset.indexOf(x) < 0;
          },
          tointeger = function(x) {
            return parseInt(x, 10);
          };
      if (control.accumulate)
        for (i=0; i <= value; i++)
          adds = adds.concat(control.subsets[i]);
      else
        adds = adds.concat(control.subsets[value]);
      deletes = fullset.filter(function(x) { return adds.indexOf(x) < 0 });  
      for (i = 0; i < subscenes.length; i++) {
        subsceneid = subscenes[i];
        if (typeof this.getObj(subsceneid) === "undefined")
          this.alertOnce("typeof object is undefined");
        for (j = 0; j < adds.length; j++)
          this.addToSubscene(adds[j], subsceneid);
        for (j = 0; j < deletes.length; j++)
          this.delFromSubscene(deletes[j], subsceneid);
      }
    };
    this.propertySetter = function(el, control)  {
      var value = control.value,
          values = [].concat(control.values),
          svals = [].concat(control.param),
          direct = values[0] === null,
          entries = [].concat(control.entries),
          ncol = entries.length,
          nrow = values.length/ncol,
          properties = this.repeatToLen(control.properties, ncol),
          objids = this.repeatToLen(control.objids, ncol),
          property, objid = objids[0],
          obj = this.getObj(objid),
          propvals, i, v1, v2, p, entry, gl, needsBinding,
          newprop, newid,
          getPropvals = function() {
            if (property === "userMatrix")
              return obj.par3d.userMatrix.getAsArray();
            else if (property === "scale" || property === "FOV" || property === "zoom")
              return [].concat(obj.par3d[property]);
            else
              return [].concat(obj[property]);
          };
          
          putPropvals = function(newvals) {
            if (newvals.length == 1)
              newvals = newvals[0];
            if (property === "userMatrix")
              obj.par3d.userMatrix.load(newvals);
            else if (property === "scale" || property === "FOV" || property === "zoom")
              obj.par3d[property] = newvals;
            else
              obj[property] = newvals;
          }
      if (direct && typeof value === "undefined")
        return;
      if (control.interp) {
        values = values.slice(0, ncol).concat(values).
                 concat(values.slice(ncol*(nrow-1), ncol*nrow));
        svals = [-Infinity].concat(svals).concat(Infinity);
        for (i = 1; i < svals.length; i++) {
          if (value <= svals[i]) {
            if (svals[i] === Infinity)
              p = 1;
            else
              p = (svals[i] - value)/(svals[i] - svals[i-1]);
            break;
          }
        }
      } else if (!direct) {
        value = Math.round(value);
      }
      for (j=0; j<entries.length; j++) {
        entry = entries[j];
        newprop = properties[j];
        newid = objids[j];
        if (newprop !== property || newid != objid) {
          if (typeof property !== "undefined")
            putPropvals(propvals);
          property = newprop;
          objid = newid;
          obj = this.getObj(objid);
          propvals = getPropvals();
        }
        if (control.interp) {
          v1 = values[ncol*(i-1) + j];
          v2 = values[ncol*i + j];
          this.setElement(propvals, entry, p*v1 + (1-p)*v2);
        } else if (!direct) {
          this.setElement(propvals, entry, values[ncol*value + j]);
        } else {
          this.setElement(propvals, entry, value[j]);
        }
      }
      putPropvals(propvals);
        
      needsBinding = [];
      for (j=0; j < entries.length; j++) {
        if (properties[j] === "values" &&
            needsBinding.indexOf(objids[j]) === -1) {
          needsBinding.push(objids[j]);
        }
      }
      for (j=0; j < needsBinding.length; j++) {
        gl = this.gl || this.initGL();
        obj = this.getObj(needsBinding[j]);
        gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
        gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
      }
    };
    this.vertexSetter = function(el, control)  {
      var svals = [].concat(control.param),
          j, k, p, propvals, stride, ofs, obj,
          attrib,
          ofss    = {x:"vofs", y:"vofs", z:"vofs",
                     red:"cofs", green:"cofs", blue:"cofs",
                     alpha:"cofs", radii:"radofs",
                     nx:"nofs", ny:"nofs", nz:"nofs",
                     ox:"oofs", oy:"oofs", oz:"oofs",
                     ts:"tofs", tt:"tofs"},
          pos     = {x:0, y:1, z:2,
                     red:0, green:1, blue:2,
                     alpha:3,radii:0,
                     nx:0, ny:1, nz:2,
                     ox:0, oy:1, oz:2,
                     ts:0, tt:1},
        values = control.values,
        direct = values === null,
        ncol,
        interp = control.interp,
        vertices = [].concat(control.vertices),
        attributes = [].concat(control.attributes),
        value = control.value;
      ncol = Math.max(vertices.length, attributes.length);
      if (!ncol)
        return;
      vertices = this.repeatToLen(vertices, ncol);
      attributes = this.repeatToLen(attributes, ncol);
      if (direct)
        interp = false;
      /* JSON doesn't pass Infinity */
      svals[0] = -Infinity;
      svals[svals.length - 1] = Infinity;
      for (j = 1; j < svals.length; j++) {
        if (value <= svals[j]) {
          if (interp) {
            if (svals[j] === Infinity)
              p = 1;
            else
              p = (svals[j] - value)/(svals[j] - svals[j-1]);
          } else {
            if (svals[j] - value > value - svals[j-1])
              j = j - 1;
          }
          break;
        }
      }
      obj = this.getObj(control.objid);
      propvals = obj.values;
      for (k=0; k<ncol; k++) {
        attrib = attributes[k];
        vertex = vertices[k];
        ofs = obj.vOffsets[ofss[attrib]];
        if (ofs < 0)
          this.alertOnce("Attribute '"+attrib+"' not found in object "+control.objid);
        else {
          stride = obj.vOffsets.stride;
          ofs = vertex*stride + ofs + pos[attrib];
          if (direct) {
            propvals[ofs] = value;
          } else if (interp) {
            propvals[ofs] = p*values[j-1][k] + (1-p)*values[j][k];
          } else {
            propvals[ofs] = values[j][k];
          }
        }
      }
      if (typeof obj.buf !== "undefined") {
        var gl = this.gl || this.initGL();
        gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
        gl.bufferData(gl.ARRAY_BUFFER, propvals, gl.STATIC_DRAW);
      }
    };
    this.ageSetter = function(el, control) {
      var objids = [].concat(control.objids),
          nobjs = objids.length,
          time = control.value,
          births = [].concat(control.births),
          ages = [].concat(control.ages),
          steps = births.length,
          j = Array(steps),
          p = Array(steps),
          i, k, age, j0, propvals, stride, ofs, objid, obj,
          attrib, dim,
          attribs = ["colors", "alpha", "radii", "vertices",
                     "normals", "origins", "texcoords",
                     "x", "y", "z",
                     "red", "green", "blue"],
          ofss    = ["cofs", "cofs", "radofs", "vofs",
                     "nofs", "oofs", "tofs",
                     "vofs", "vofs", "vofs",
                     "cofs", "cofs", "cofs"],
          dims    = [3,1,1,3,
                     3,2,2,
                     1,1,1,
                     1,1,1],
          pos     = [0,3,0,0,
                     0,0,0,
                     0,1,2,
                     0,1,2];
      /* Infinity doesn't make it through JSON */
      ages[0] = -Infinity;
      ages[ages.length-1] = Infinity;
      for (i = 0; i < steps; i++) {
        if (births[i] !== null) {  // NA in R becomes null
          age = time - births[i];
          for (j0 = 1; age > ages[j0]; j0++);
          if (ages[j0] == Infinity)
            p[i] = 1;
          else if (ages[j0] > ages[j0-1])
            p[i] = (ages[j0] - age)/(ages[j0] - ages[j0-1]);
          else
            p[i] = 0;
          j[i] = j0;
        }
      }
      for (l = 0; l < nobjs; l++) {
        objid = objids[l];
        obj = this.getObj(objid);
        if (typeof obj.vOffsets === "undefined")
          continue;
        propvals = obj.values;
        stride = obj.vOffsets.stride;
        for (k = 0; k < attribs.length; k++) {
          attrib = control[attribs[k]];
          if (typeof attrib !== "undefined") {
            ofs = obj.vOffsets[ofss[k]];
            if (ofs >= 0) {
              dim = dims[k];
              ofs = ofs + pos[k];
              for (i = 0; i < steps; i++) {
                if (births[i] !== null) {
                  for (d=0; d < dim; d++) {
                    propvals[i*stride + ofs + d] = p[i]*attrib[dim*(j[i]-1) + d] + (1-p[i])*attrib[dim*j[i] + d];
                  }
                }
              }
            } else
              this.alertOnce("\'"+attribs[k]+"\' property not found in object "+objid);
          }
        }
        obj.values = propvals;
        if (typeof obj.buf !== "undefined") {
          gl = this.gl || this.initGL();
          gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
          gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
        }
      }
    };
    this.oldBridge = function(el, control) {
      var attrname, global = window[control.prefix + "rgl"];
      if (typeof global !== "undefined")
        for (attrname in global)
          this[attrname] = global[attrname];
      window[control.prefix + "rgl"] = this;
    };
    this.Player = function(el, control) {
      var
        self = this,
        components = [].concat(control.components),
        buttonLabels = [].concat(control.buttonLabels),
        Tick = function() { /* "this" will be a timer */
          var i,
              nominal = this.value,
              slider = this.Slider,
              labels = this.outputLabels,
              output = this.Output,
              step;
          if (typeof slider !== "undefined" && nominal != slider.value)
            slider.value = nominal;
          if (typeof output !== "undefined") {
            step = Math.round((nominal - output.sliderMin)/output.sliderStep);
            if (labels !== null) {
              output.innerHTML = labels[step];
            } else {
              step = step*output.sliderStep + output.sliderMin;
              output.innerHTML = step.toPrecision(output.outputPrecision);
            }
          }
          for (i=0; i < this.actions.length; i++) {
            this.actions[i].value = nominal;
          }
          self.applyControls(el, this.actions, false);
          self.drawScene();
        },
        OnSliderInput = function() { /* "this" will be the slider */
          this.rgltimer.value = Number(this.value);
          this.rgltimer.Tick();
        },
        addSlider = function(min, max, step, value) {
          var slider = document.createElement("input");
          slider.type = "range";
          slider.min = min;
          slider.max = max;
          slider.step = step;
          slider.value = value;
          slider.oninput = OnSliderInput;
          slider.sliderActions = control.actions;
          slider.sliderScene = this;
          slider.className = "rgl-slider";
          slider.id = el.id + "-slider";
          el.rgltimer.Slider = slider;
          slider.rgltimer = el.rgltimer;
          el.appendChild(slider);
        },
        addLabel = function(labels, min, step, precision) {
          var output = document.createElement("output");
          output.sliderMin = min;
          output.sliderStep = step;
          output.outputPrecision = precision;
          output.className = "rgl-label";
          output.id = el.id + "-label";
          el.rgltimer.Output = output;
          el.rgltimer.outputLabels = labels;
          el.appendChild(output);
        },
        addButton = function(which, label, active) {
          var button = document.createElement("input"),
              onclicks = {Reverse: function() { this.rgltimer.reverse();},
                    Play: function() { this.rgltimer.play();
                                       this.value = this.rgltimer.enabled ? this.inactiveValue : this.activeValue; },
                   Slower: function() { this.rgltimer.slower(); },
                   Faster: function() { this.rgltimer.faster(); },
                   Reset: function() { this.rgltimer.reset(); },
              	   Step:  function() { this.rgltimer.step(); }
              };
          button.rgltimer = el.rgltimer;
          button.type = "button";
          button.value = label;
          button.activeValue = label;
          button.inactiveValue = active;
          if (which === "Play")
            button.rgltimer.PlayButton = button;
          button.onclick = onclicks[which];
          button.className = "rgl-button";
          button.id = el.id + "-" + which;
          el.appendChild(button);
        };
        if (typeof control.reinit !== "undefined" && control.reinit !== null) {
          control.actions.reinit = control.reinit;
        }
        el.rgltimer = new rgltimerClass(Tick, control.start, control.interval, control.stop, 
                                        control.step, control.value, control.rate, control.loop, control.actions);
        for (var i=0; i < components.length; i++) {
          switch(components[i]) {
            case "Slider": addSlider(control.start, control.stop,
                                   control.step, control.value);
              break;
            case "Label": addLabel(control.labels, control.start,
                                   control.step, control.precision);
              break;
            default:
              addButton(components[i], buttonLabels[i], control.pause);
          }
        }
        el.rgltimer.Tick();
    };
    this.applyControls = function(el, x, draw) {
      var self = this, reinit = x.reinit, i, control, type;
      for (i = 0; i < x.length; i++) {
        control = x[i];
        type = control.type;
        self[type](el, control);
      }
      if (typeof reinit !== "undefined" && reinit !== null) {
        reinit = [].concat(reinit);
        for (i = 0; i < reinit.length; i++)
          self.getObj(reinit[i]).initialized = false;
      }
      if (typeof draw === "undefined" || draw)
        self.drawScene();
    };
    this.sceneChangeHandler = function(message) {
      var self = document.getElementById(message.elementId).rglinstance,
          objs = message.objects, mat = message.material,
          root = message.rootSubscene,
          initSubs = message.initSubscenes,
          redraw = message.redrawScene,
          skipRedraw = message.skipRedraw,
          deletes, subs, allsubs = [], i,j;
      if (typeof message.delete !== "undefined") {
        deletes = [].concat(message.delete);
        if (typeof message.delfromSubscenes !== "undefined")
          subs = [].concat(message.delfromSubscenes);
        else
          subs = [];
        for (i = 0; i < deletes.length; i++) {
          for (j = 0; j < subs.length; j++) {
            self.delFromSubscene(deletes[i], subs[j]);
          }
          delete self.scene.objects[deletes[i]];
        }
      }
      if (typeof objs !== "undefined") {
        Object.keys(objs).forEach(function(key){
          key = parseInt(key, 10);
          self.scene.objects[key] = objs[key];
          self.initObj(key);
          var obj = self.getObj(key),
              subs = [].concat(obj.inSubscenes), k;
          allsubs = allsubs.concat(subs);
          for (k = 0; k < subs.length; k++)
            self.addToSubscene(key, subs[k]);
        });
      }
      if (typeof mat !== "undefined") {
        self.scene.material = mat;
      }
      if (typeof root !== "undefined") {
        self.scene.rootSubscene = root;
      }
      if (typeof initSubs !== "undefined")
        allsubs = allsubs.concat(initSubs);
      allsubs = self.unique(allsubs);
      for (i = 0; i < allsubs.length; i++) {
        self.initSubscene(allsubs[i]);
      }
      if (typeof skipRedraw !== "undefined") {
        root = self.getObj(self.scene.rootSubscene);
        root.par3d.skipRedraw = skipRedraw;
      }
      if (redraw)
        self.drawScene();
    };
}).call(rglwidgetClass.prototype);
rgltimerClass = function(Tick, startTime, interval, stopTime, stepSize, value, rate, loop, actions) {
  this.enabled = false;
  this.timerId = 0;
  this.startTime = startTime;         /* nominal start time in seconds */
  this.value = value;                 /* current nominal time */
  this.interval = interval;           /* seconds between updates */
  this.stopTime = stopTime;           /* nominal stop time */
  this.stepSize = stepSize;           /* nominal step size */
  this.rate = rate;                   /* nominal units per second */
  this.loop = loop;                   /* "none", "cycle", or "oscillate" */
  this.realStart = undefined;         /* real world start time */
  this.multiplier = 1;                /* multiplier for fast-forward
                                         or reverse */
  this.actions = actions;
  this.Tick = Tick;
};
(function() {
  this.play = function() {
    if (this.enabled) {
      this.enabled = false;
      window.clearInterval(this.timerId);
      this.timerId = 0;
      return;
    }
    var tick = function(self) {
      var now = new Date();
      self.value = self.multiplier*self.rate*(now - self.realStart)/1000 + self.startTime;
      self.forceToRange();
      if (typeof self.Tick !== "undefined") {
        self.Tick(self.value);
      }
    };
    this.realStart = new Date() - 1000*(this.value - this.startTime)/this.rate/this.multiplier;
    this.timerId = window.setInterval(tick, 1000*this.interval, this);
    this.enabled = true;
  };
  
  this.forceToRange = function() {
    if (this.value > this.stopTime + this.stepSize/2 || this.value < this.startTime - this.stepSize/2) {
      if (!this.loop) {
        this.reset();
      } else {
        var cycle = this.stopTime - this.startTime + this.stepSize,
            newval = (this.value - this.startTime) % cycle + this.startTime;
        if (newval < this.startTime) {
          newval += cycle;
        }
        this.realStart += (this.value - newval)*1000/this.multiplier/this.rate;
        this.value = newval;
      }
    }  	
  }
  this.reset = function() {
    this.value = this.startTime;
    this.newmultiplier(1);
    if (typeof this.Tick !== "undefined") {
        this.Tick(this.value);
    }
    if (this.enabled)
      this.play();  /* really pause... */
    if (typeof this.PlayButton !== "undefined")
      this.PlayButton.value = "Play";
  };
  this.faster = function() {
    this.newmultiplier(Math.SQRT2*this.multiplier);
  };
  this.slower = function() {
    this.newmultiplier(this.multiplier/Math.SQRT2);
  };
  this.reverse = function() {
    this.newmultiplier(-this.multiplier);
  };
  this.newmultiplier = function(newmult) {
    if (newmult != this.multiplier) {
      this.realStart += 1000*(this.value - this.startTime)/this.rate*(1/this.multiplier - 1/newmult);
      this.multiplier = newmult;
    }
  };
  
  this.step = function() {
    this.value += this.rate*this.multiplier;
    this.forceToRange();
    if (typeof this.Tick !== "undefined")
      this.Tick(this.value);
  }
}).call(rgltimerClass.prototype);</script>

<div id="div" class="rglWebGL"></div>
<script type="text/javascript">
	var div = document.getElementById("div"),
      rgl = new rglwidgetClass();
  div.width = 1400;
  div.height = 719.1257;
  rgl.initialize(div,
                         {"material":{"color":"#000000","alpha":1,"lit":true,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":false,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less","isTransparent":false},"rootSubscene":115,"objects":{"169":{"id":169,"type":"linestrip","material":{"lit":false},"vertices":[[-5,-5,1.200921e-005],[-4.5,-5,1.817341e-005],[-4,-5,2.672699e-005],[-3.5,-5,3.819928e-005],[-3,-5,5.305813e-005],[-2.5,-5,7.162098e-005],[-2,-5,9.395506e-005],[-1.5,-5,0.000119782],[-1,-5,0.0001484071],[-0.5,-5,0.0001786938],[0,-5,0.0002091008],[0.5,-5,0.00023779],[1,-5,0.0002627986],[1.5,-5,0.0002822566],[2,-5,0.0002946163],[2.5,-5,0.0002988553],[3,-5,0.0002946163],[3.5,-5,0.0002822566],[4,-5,0.0002627986],[4.5,-5,0.00023779],[5,-5,0.0002091008],[5.5,-5,0.0001786938],[6,-5,0.0001484071],[6.5,-5,0.000119782],[7,-5,9.395506e-005],[7.5,-5,7.162098e-005],[8,-5,5.305813e-005],[8.5,-5,3.819928e-005],[9,-5,2.672699e-005],[9.5,-5,1.817341e-005],[10,-5,1.200921e-005],[10.5,-5,7.712301e-006],[11,-5,4.813325e-006],[11.5,-5,2.919429e-006],[12,-5,1.720847e-006],[12.5,-5,9.857758e-007],[13,-5,5.487894e-007],[13.5,-5,2.9691e-007],[14,-5,1.561117e-007],[14.5,-5,7.976966e-008],[15,-5,3.961244e-008],[-5,-4.5,2.211792e-005],[-4.5,-4.5,3.371074e-005],[-4,-4.5,4.993258e-005],[-3.5,-4.5,7.187723e-005],[-3,-4.5,0.0001005519],[-2.5,-4.5,0.0001367038],[-2,-4.5,0.0001806187],[-1.5,-4.5,0.000231919],[-1,-4.5,0.000289402],[-0.5,-4.5,0.0003509605],[0,-4.5,0.0004136249],[0.5,-4.5,0.0004737471],[1,-4.5,0.0005273247],[1.5,-4.5,0.0005704286],[2,-4.5,0.0005996751],[2.5,-4.5,0.000612664],[3,-4.5,0.0006083034],[3.5,-4.5,0.0005869616],[4,-4.5,0.0005504157],[4.5,-4.5,0.0005016068],[5,-4.5,0.0004442502],[5.5,-4.5,0.0003823697],[6,-4.5,0.0003198386],[6.5,-4.5,0.000259998],[7,-4.5,0.0002054],[7.5,-4.5,0.0001576967],[8,-4.5,0.0001176621],[8.5,-4.5,8.531822e-005],[9,-4.5,6.012273e-005],[9.5,-4.5,4.11744e-005],[10,-4.5,2.740358e-005],[10.5,-4.5,1.772471e-005],[11,-4.5,1.114146e-005],[11.5,-4.5,6.806082e-006],[12,-4.5,4.040579e-006],[12.5,-4.5,2.331211e-006],[13,-4.5,1.307107e-006],[13.5,-4.5,7.1225e-007],[14,-4.5,3.771771e-007],[14.5,-4.5,1.941108e-007],[15,-4.5,9.708354e-008],[-5,-4,3.819928e-005],[-4.5,-4,5.86383e-005],[-4,-4,8.747806e-005],[-3.5,-4,0.0001268261],[-3,-4,0.0001786938],[-2.5,-4,0.000244682],[-2,-4,0.0003256013],[-1.5,-4,0.0004210774],[-1,-4,0.0005292114],[-0.5,-4,0.0006463803],[0,-4,0.000767253],[0.5,-4,0.0008850762],[1,-4,0.0009922344],[1.5,-4,0.001081035],[2,-4,0.001144607],[2.5,-4,0.001177782],[3,-4,0.001177782],[3.5,-4,0.001144607],[4,-4,0.001081035],[4.5,-4,0.0009922344],[5,-4,0.0008850762],[5.5,-4,0.000767253],[6,-4,0.0006463803],[6.5,-4,0.0005292114],[7,-4,0.0004210774],[7.5,-4,0.0003256013],[8,-4,0.000244682],[8.5,-4,0.0001786938],[9,-4,0.0001268261],[9.5,-4,8.747806e-005],[10,-4,5.86383e-005],[10.5,-4,3.819928e-005],[11,-4,2.418358e-005],[11.5,-4,1.487913e-005],[12,-4,8.896641e-006],[12.5,-4,5.16971e-006],[13,-4,2.919429e-006],[13.5,-4,1.602217e-006],[14,-4,8.545474e-007],[14.5,-4,4.429376e-007],[15,-4,2.231211e-007],[-5,-3.5,6.18653e-005],[-4.5,-3.5,9.56479e-005],[-4,-3.5,0.0001437128],[-3.5,-3.5,0.000209849],[-3,-3.5,0.0002977898],[-2.5,-3.5,0.000410681],[-2,-3.5,0.0005504157],[-1.5,-3.5,0.0007169166],[-1,-3.5,0.000907482],[-0.5,-3.5,0.001116347],[0,-3.5,0.001334602],[0.5,-3.5,0.001550586],[1,-3.5,0.00175078],[1.5,-3.5,0.001921139],[2,-3.5,0.002048697],[2.5,-3.5,0.002123188],[3,-3.5,0.002138407],[3.5,-3.5,0.002093072],[4,-3.5,0.001990991],[4.5,-3.5,0.001840544],[5,-3.5,0.00165354],[5.5,-3.5,0.001443693],[6,-3.5,0.001224973],[6.5,-3.5,0.001010112],[7,-3.5,0.0008094767],[7.5,-3.5,0.0006304211],[8,-3.5,0.0004771432],[8.5,-3.5,0.0003509605],[9,-3.5,0.0002508762],[9.5,-3.5,0.0001742818],[10,-3.5,0.0001176621],[10.5,-3.5,7.719912e-005],[11,-3.5,4.922433e-005],[11.5,-3.5,3.050274e-005],[12,-3.5,1.836917e-005],[12.5,-3.5,1.075058e-005],[13,-3.5,6.114566e-006],[13.5,-3.5,3.3798e-006],[14,-3.5,1.815549e-006],[14.5,-3.5,9.477998e-007],[15,-3.5,4.808579e-007],[-5,-3,9.395506e-005],[-4.5,-3,0.0001463021],[-4,-3,0.0002213974],[-3.5,-3,0.0003256013],[-3,-3,0.0004653625],[-2.5,-3,0.0006463803],[-2,-3,0.0008725221],[-1.5,-3,0.001144607],[-1,-3,0.001459244],[-0.5,-3,0.001807969],[0,-3,0.002176936],[0.5,-3,0.00254737],[1,-3,0.002896876],[1.5,-3,0.003201543],[2,-3,0.003438589],[2.5,-3,0.003589161],[3,-3,0.003640803],[3.5,-3,0.003589161],[4,-3,0.003438589],[4.5,-3,0.003201543],[5,-3,0.002896876],[5.5,-3,0.00254737],[6,-3,0.002176936],[6.5,-3,0.001807969],[7,-3,0.001459244],[7.5,-3,0.001144607],[8,-3,0.0008725221],[8.5,-3,0.0006463803],[9,-3,0.0004653625],[9.5,-3,0.0003256013],[10,-3,0.0002213974],[10.5,-3,0.0001463021],[11,-3,9.395506e-005],[11.5,-3,5.86383e-005],[12,-3,3.556592e-005],[12.5,-3,2.096421e-005],[13,-3,1.200921e-005],[13.5,-3,6.685623e-006],[14,-3,3.617104e-006],[14.5,-3,1.90183e-006],[15,-3,9.717935e-007],[-5,-2.5,0.0001338056],[-4.5,-2.5,0.000209849],[-4,-2.5,0.0003198386],[-3.5,-2.5,0.0004737471],[-3,-2.5,0.0006819521],[-2.5,-2.5,0.0009540096],[-2,-2.5,0.00129701],[-1.5,-2.5,0.001713662],[-1,-2.5,0.002200386],[-0.5,-2.5,0.002745769],[0,-2.5,0.003329821],[0.5,-2.5,0.003924364],[1,-2.5,0.004494788],[1.5,-2.5,0.005003118],[2,-2.5,0.005412076],[2.5,-2.5,0.005689559],[3,-2.5,0.005812794],[3.5,-2.5,0.005771422],[4,-2.5,0.005568937],[4.5,-2.5,0.005222199],[5,-2.5,0.004759113],[5.5,-2.5,0.004214929],[6,-2.5,0.003627823],[6.5,-2.5,0.003034544],[7,-2.5,0.002466792],[7.5,-2.5,0.001948781],[8,-2.5,0.001496185],[8.5,-2.5,0.001116347],[9,-2.5,0.0008094767],[9.5,-2.5,0.0005704286],[10,-2.5,0.0003906518],[10.5,-2.5,0.000259998],[11,-2.5,0.0001681673],[11.5,-2.5,0.0001057073],[12,-2.5,6.457431e-005],[12.5,-2.5,3.833595e-005],[13,-2.5,2.211792e-005],[13.5,-2.5,1.240149e-005],[14,-2.5,6.75764e-006],[14.5,-2.5,3.578556e-006],[15,-2.5,1.841672e-006],[-5,-2,0.0001786938],[-4.5,-2,0.0002822566],[-4,-2,0.0004332817],[-3.5,-2,0.0006463803],[-3,-2,0.0009371249],[-2.5,-2,0.001320379],[-2,-2,0.001807969],[-1.5,-2,0.002405887],[-1,-2,0.003111364],[-0.5,-2,0.003910373],[0,-2,0.00477614],[0.5,-2,0.005669275],[1,-2,0.006539878],[1.5,-2,0.007331676],[2,-2,0.007987824],[2.5,-2,0.008457565],[3,-2,0.008702695],[3.5,-2,0.008702695],[4,-2,0.008457565],[4.5,-2,0.007987824],[5,-2,0.007331676],[5.5,-2,0.006539878],[6,-2,0.005669275],[6.5,-2,0.00477614],[7,-2,0.003910373],[7.5,-2,0.003111364],[8,-2,0.002405887],[8.5,-2,0.001807969],[9,-2,0.001320379],[9.5,-2,0.0009371249],[10,-2,0.0006463803],[10.5,-2,0.0004332817],[11,-2,0.0002822566],[11.5,-2,0.0001786938],[12,-2,0.0001099427],[12.5,-2,6.573778e-005],[13,-2,3.819928e-005],[13.5,-2,2.157182e-005],[14,-2,1.183887e-005],[14.5,-2,6.314298e-006],[15,-2,3.272891e-006],[-5,-1.5,0.0002237823],[-4.5,-1.5,0.0003560102],[-4,-1.5,0.0005504157],[-3.5,-1.5,0.0008270098],[-3,-1.5,0.001207597],[-2.5,-1.5,0.001713662],[-2,-1.5,0.002363306],[-1.5,-1.5,0.003167423],[-1,-1.5,0.00412557],[-0.5,-1.5,0.005222199],[0,-1.5,0.006424131],[0.5,-1.5,0.007680101],[1,-1.5,0.008923004],[1.5,-1.5,0.01007504],[2,-1.5,0.01105539],[2.5,-1.5,0.01178944],[3,-1.5,0.0122181],[3.5,-1.5,0.01230568],[4,-1.5,0.0120448],[4.5,-1.5,0.01145737],[5,-1.5,0.0105916],[5.5,-1.5,0.009515466],[6,-1.5,0.008307878],[6.5,-1.5,0.00704923],[7,-1.5,0.005812794],[7.5,-1.5,0.004658217],[8,-1.5,0.003627823],[8.5,-1.5,0.002745769],[9,-1.5,0.002019638],[9.5,-1.5,0.001443693],[10,-1.5,0.001002923],[10.5,-1.5,0.0006770984],[11,-1.5,0.0004442502],[11.5,-1.5,0.0002832665],[12,-1.5,0.0001755312],[12.5,-1.5,0.0001057073],[13,-1.5,6.18653e-005],[13.5,-1.5,3.518689e-005],[14,-1.5,1.944941e-005],[14.5,-1.5,1.044776e-005],[15,-1.5,5.454212e-006],[-5,-1,0.0002627986],[-4.5,-1,0.0004210774],[-4,-1,0.0006556806],[-3.5,-1,0.0009922344],[-3,-1,0.001459244],[-2.5,-1,0.00208561],[-2,-1,0.002896876],[-1.5,-1,0.003910373],[-1,-1,0.005129773],[-0.5,-1,0.006539878],[0,-1,0.008102755],[0.5,-1,0.009756351],[1,-1,0.01141652],[1.5,-1,0.0129829],[2,-1,0.01434832],[2.5,-1,0.01541069],[3,-1,0.0160855],[3.5,-1,0.01631695],[4,-1,0.0160855],[4.5,-1,0.01541069],[5,-1,0.01434832],[5.5,-1,0.0129829],[6,-1,0.01141652],[6.5,-1,0.009756351],[7,-1,0.008102755],[7.5,-1,0.006539878],[8,-1,0.005129773],[8.5,-1,0.003910373],[9,-1,0.002896876],[9.5,-1,0.00208561],[10,-1,0.001459244],[10.5,-1,0.0009922344],[11,-1,0.0006556806],[11.5,-1,0.0004210774],[12,-1,0.0002627986],[12.5,-1,0.0001593954],[13,-1,9.395506e-005],[13.5,-1,5.382154e-005],[14,-1,2.996288e-005],[14.5,-1,1.621074e-005],[15,-1,8.523411e-006],[-5,-0.5,0.000289402],[-4.5,-0.5,0.0004670275],[-4,-0.5,0.0007324449],[-3.5,-0.5,0.001116347],[-3,-0.5,0.00165354],[-2.5,-0.5,0.002380247],[-2,-0.5,0.003329821],[-1.5,-0.5,0.004527009],[-1,-0.5,0.005981269],[-0.5,-0.5,0.007680101],[0,-0.5,0.009583677],[0.5,-0.5,0.01162222],[1,-0.5,0.01369737],[1.5,-0.5,0.01568835],[2,-0.5,0.0174626],[2.5,-0.5,0.01889],[3,-0.5,0.01985851],[3.5,-0.5,0.02028864],[4,-0.5,0.02014424],[4.5,-0.5,0.0194375],[5,-0.5,0.01822726],[5.5,-0.5,0.01661094],[6,-0.5,0.01471155],[6.5,-0.5,0.01266235],[7,-0.5,0.0105916],[7.5,-0.5,0.00860995],[8,-0.5,0.006801915],[8.5,-0.5,0.005222199],[9,-0.5,0.003896432],[9.5,-0.5,0.002825351],[10,-0.5,0.001990991],[10.5,-0.5,0.001363509],[11,-0.5,0.000907482],[11.5,-0.5,0.0005869616],[12,-0.5,0.0003689546],[12.5,-0.5,0.0002253865],[13,-0.5,0.0001338056],[13.5,-0.5,7.719912e-005],[14,-0.5,4.328545e-005],[14.5,-0.5,2.358648e-005],[15,-0.5,1.249039e-005],[-5,0,0.0002988553],[-4.5,0,0.0004857401],[-4,0,0.000767253],[-3.5,0,0.001177782],[-3,0,0.001757044],[-2.5,0,0.00254737],[-2,0,0.003589161],[-1.5,0,0.00491457],[-1,0,0.006539878],[-0.5,0,0.008457565],[0,0,0.0106295],[0.5,0,0.0129829],[1,0,0.01541069],[1.5,0,0.01777723],[2,0,0.01992956],[2.5,0,0.02171316],[3,0,0.02299005],[3.5,0,0.02365638],[4,0,0.02365638],[4.5,0,0.02299005],[5,0,0.02171316],[5.5,0,0.01992956],[6,0,0.01777723],[6.5,0,0.01541069],[7,0,0.0129829],[7.5,0,0.0106295],[8,0,0.008457565],[8.5,0,0.006539878],[9,0,0.00491457],[9.5,0,0.003589161],[10,0,0.00254737],[10.5,0,0.001757044],[11,0,0.001177782],[11.5,0,0.000767253],[12,0,0.0004857401],[12.5,0,0.0002988553],[13,0,0.0001786938],[13.5,0,0.0001038364],[14,0,5.86383e-005],[14.5,0,3.218138e-005],[15,0,1.716404e-005],[-5,0.5,0.000289402],[-4.5,0.5,0.0004737471],[-4,0.5,0.0007536737],[-3.5,0.5,0.00116523],[-3,0.5,0.00175078],[-2.5,0.5,0.002556484],[-2,0.5,0.003627823],[-1.5,0.5,0.005003118],[-1,0.5,0.006705435],[-0.5,0.5,0.008733831],[0,0.5,0.01105539],[0.5,0.5,0.01359988],[1,0.5,0.01625877],[1.5,0.5,0.01889],[2,0.5,0.02132886],[2.5,0.5,0.02340427],[3,0.5,0.02495824],[3.5,0.5,0.02586572],[4,0.5,0.02605113],[4.5,0.5,0.02549884],[5,0.5,0.02425524],[5.5,0.5,0.02242242],[6,0.5,0.02014424],[6.5,0.5,0.01758778],[7,0.5,0.01492322],[7.5,0.5,0.01230568],[8,0.5,0.009861445],[8.5,0.5,0.007680101],[9,0.5,0.005812794],[9.5,0.5,0.004275574],[10,0.5,0.003056297],[10.5,0.5,0.002123188],[11,0.5,0.001433417],[11.5,0.5,0.0009404778],[12,0.5,0.0005996751],[12.5,0.5,0.0003715995],[13,0.5,0.0002237823],[13.5,0.5,0.0001309688],[14,0.5,7.449066e-005],[14.5,0.5,4.11744e-005],[15,0.5,2.211792e-005],[-5,1,0.0002627986],[-4.5,1,0.0004332817],[-4,1,0.0006942392],[-3.5,1,0.001081035],[-3,1,0.001635918],[-2.5,1,0.002405887],[-2,1,0.003438589],[-1.5,1,0.00477614],[-1,1,0.006447115],[-0.5,1,0.008457565],[0,1,0.01078244],[0.5,1,0.01335919],[1,1,0.0160855],[1.5,1,0.01882266],[2,1,0.02140518],[2.5,1,0.02365638],[3,1,0.02540793],[3.5,1,0.02652051],[4,1,0.0269021],[4.5,1,0.02652051],[5,1,0.02540793],[5.5,1,0.02365638],[6,1,0.02140518],[6.5,1,0.01882266],[7,1,0.0160855],[7.5,1,0.01335919],[8,1,0.01078244],[8.5,1,0.008457565],[9,1,0.006447115],[9.5,1,0.00477614],[10,1,0.003438589],[10.5,1,0.002405887],[11,1,0.001635918],[11.5,1,0.001081035],[12,1,0.0006942392],[12.5,1,0.0004332817],[13,1,0.0002627986],[13.5,1,0.0001549057],[14,1,8.873671e-005],[14.5,1,4.940045e-005],[15,1,2.672699e-005],[-5,1.5,0.0002237823],[-4.5,1.5,0.0003715995],[-4,1.5,0.0005996751],[-3.5,1.5,0.0009404778],[-3,1.5,0.001433417],[-2.5,1.5,0.002123188],[-2,1.5,0.003056297],[-1.5,1.5,0.004275574],[-1,1.5,0.005812794],[-0.5,1.5,0.007680101],[0,1.5,0.009861445],[0.5,1.5,0.01230568],[1,1.5,0.01492322],[1.5,1.5,0.01758778],[2,1.5,0.02014424],[2.5,1.5,0.02242242],[3,1.5,0.02425524],[3.5,1.5,0.02549884],[4,1.5,0.02605113],[4.5,1.5,0.02586572],[5,1.5,0.02495824],[5.5,1.5,0.02340427],[6,1.5,0.02132886],[6.5,1.5,0.01889],[7,1.5,0.01625877],[7.5,1.5,0.01359988],[8,1.5,0.01105539],[8.5,1.5,0.008733831],[9,1.5,0.006705435],[9.5,1.5,0.005003118],[10,1.5,0.003627823],[10.5,1.5,0.002556484],[11,1.5,0.00175078],[11.5,1.5,0.00116523],[12,1.5,0.0007536737],[12.5,1.5,0.0004737471],[13,1.5,0.000289402],[13.5,1.5,0.0001718098],[14,1.5,9.912563e-005],[14.5,1.5,5.557962e-005],[15,1.5,3.028564e-005],[-5,2,0.0001786938],[-4.5,2,0.0002988553],[-4,2,0.0004857401],[-3.5,2,0.000767253],[-3,2,0.001177782],[-2.5,2,0.001757044],[-2,2,0.00254737],[-1.5,2,0.003589161],[-1,2,0.00491457],[-0.5,2,0.006539878],[0,2,0.008457565],[0.5,2,0.0106295],[1,2,0.0129829],[1.5,2,0.01541069],[2,2,0.01777723],[2.5,2,0.01992956],[3,2,0.02171316],[3.5,2,0.02299005],[4,2,0.02365638],[4.5,2,0.02365638],[5,2,0.02299005],[5.5,2,0.02171316],[6,2,0.01992956],[6.5,2,0.01777723],[7,2,0.01541069],[7.5,2,0.0129829],[8,2,0.0106295],[8.5,2,0.008457565],[9,2,0.006539878],[9.5,2,0.00491457],[10,2,0.003589161],[10.5,2,0.00254737],[11,2,0.001757044],[11.5,2,0.001177782],[12,2,0.000767253],[12.5,2,0.0004857401],[13,2,0.0002988553],[13.5,2,0.0001786938],[14,2,0.0001038364],[14.5,2,5.86383e-005],[15,2,3.218138e-005],[-5,2.5,0.0001338056],[-4.5,2.5,0.0002253865],[-4,2.5,0.0003689546],[-3.5,2.5,0.0005869616],[-3,2.5,0.000907482],[-2.5,2.5,0.001363509],[-2,2.5,0.001990991],[-1.5,2.5,0.002825351],[-1,2.5,0.003896432],[-0.5,2.5,0.005222199],[0,2.5,0.006801915],[0.5,2.5,0.00860995],[1,2.5,0.0105916],[1.5,2.5,0.01266235],[2,2.5,0.01471155],[2.5,2.5,0.01661094],[3,2.5,0.01822726],[3.5,2.5,0.0194375],[4,2.5,0.02014424],[4.5,2.5,0.02028864],[5,2.5,0.01985851],[5.5,2.5,0.01889],[6,2.5,0.0174626],[6.5,2.5,0.01568835],[7,2.5,0.01369737],[7.5,2.5,0.01162222],[8,2.5,0.009583677],[8.5,2.5,0.007680101],[9,2.5,0.005981269],[9.5,2.5,0.004527009],[10,2.5,0.003329821],[10.5,2.5,0.002380247],[11,2.5,0.00165354],[11.5,2.5,0.001116347],[12,2.5,0.0007324449],[12.5,2.5,0.0004670275],[13,2.5,0.000289402],[13.5,2.5,0.0001742818],[14,2.5,0.0001019986],[14.5,2.5,5.801339e-005],[15,2.5,3.206665e-005],[-5,3,9.395506e-005],[-4.5,3,0.0001593954],[-4,3,0.0002627986],[-3.5,3,0.0004210774],[-3,3,0.0006556806],[-2.5,3,0.0009922344],[-2,3,0.001459244],[-1.5,3,0.00208561],[-1,3,0.002896876],[-0.5,3,0.003910373],[0,3,0.005129773],[0.5,3,0.006539878],[1,3,0.008102755],[1.5,3,0.009756351],[2,3,0.01141652],[2.5,3,0.0129829],[3,3,0.01434832],[3.5,3,0.01541069],[4,3,0.0160855],[4.5,3,0.01631695],[5,3,0.0160855],[5.5,3,0.01541069],[6,3,0.01434832],[6.5,3,0.0129829],[7,3,0.01141652],[7.5,3,0.009756351],[8,3,0.008102755],[8.5,3,0.006539878],[9,3,0.005129773],[9.5,3,0.003910373],[10,3,0.002896876],[10.5,3,0.00208561],[11,3,0.001459244],[11.5,3,0.0009922344],[12,3,0.0006556806],[12.5,3,0.0004210774],[13,3,0.0002627986],[13.5,3,0.0001593954],[14,3,9.395506e-005],[14.5,3,5.382154e-005],[15,3,2.996288e-005],[-5,3.5,6.18653e-005],[-4.5,3.5,0.0001057073],[-4,3.5,0.0001755312],[-3.5,3.5,0.0002832665],[-3,3.5,0.0004442502],[-2.5,3.5,0.0006770984],[-2,3.5,0.001002923],[-1.5,3.5,0.001443693],[-1,3.5,0.002019638],[-0.5,3.5,0.002745769],[0,3.5,0.003627823],[0.5,3.5,0.004658217],[1,3.5,0.005812794],[1.5,3.5,0.00704923],[2,3.5,0.008307878],[2.5,3.5,0.009515466],[3,3.5,0.0105916],[3.5,3.5,0.01145737],[4,3.5,0.0120448],[4.5,3.5,0.01230568],[5,3.5,0.0122181],[5.5,3.5,0.01178944],[6,3.5,0.01105539],[6.5,3.5,0.01007504],[7,3.5,0.008923004],[7.5,3.5,0.007680101],[8,3.5,0.006424131],[8.5,3.5,0.005222199],[9,3.5,0.00412557],[9.5,3.5,0.003167423],[10,3.5,0.002363306],[10.5,3.5,0.001713662],[11,3.5,0.001207597],[11.5,3.5,0.0008270098],[12,3.5,0.0005504157],[12.5,3.5,0.0003560102],[13,3.5,0.0002237823],[13.5,3.5,0.0001367038],[14,3.5,8.11572e-005],[14.5,3.5,4.682363e-005],[15,3.5,2.625395e-005],[-5,4,3.819928e-005],[-4.5,4,6.573778e-005],[-4,4,0.0001099427],[-3.5,4,0.0001786938],[-3,4,0.0002822566],[-2.5,4,0.0004332817],[-2,4,0.0006463803],[-1.5,4,0.0009371249],[-1,4,0.001320379],[-0.5,4,0.001807969],[0,4,0.002405887],[0.5,4,0.003111364],[1,4,0.003910373],[1.5,4,0.00477614],[2,4,0.005669275],[2.5,4,0.006539878],[3,4,0.007331676],[3.5,4,0.007987824],[4,4,0.008457565],[4.5,4,0.008702695],[5,4,0.008702695],[5.5,4,0.008457565],[6,4,0.007987824],[6.5,4,0.007331676],[7,4,0.006539878],[7.5,4,0.005669275],[8,4,0.00477614],[8.5,4,0.003910373],[9,4,0.003111364],[9.5,4,0.002405887],[10,4,0.001807969],[10.5,4,0.001320379],[11,4,0.0009371249],[11.5,4,0.0006463803],[12,4,0.0004332817],[12.5,4,0.0002822566],[13,4,0.0001786938],[13.5,4,0.0001099427],[14,4,6.573778e-005],[14.5,4,3.819928e-005],[15,4,2.157182e-005],[-5,4.5,2.211792e-005],[-4.5,4.5,3.833595e-005],[-4,4.5,6.457431e-005],[-3.5,4.5,0.0001057073],[-3,4.5,0.0001681673],[-2.5,4.5,0.000259998],[-2,4.5,0.0003906518],[-1.5,4.5,0.0005704286],[-1,4.5,0.0008094767],[-0.5,4.5,0.001116347],[0,4.5,0.001496185],[0.5,4.5,0.001948781],[1,4.5,0.002466792],[1.5,4.5,0.003034544],[2,4.5,0.003627823],[2.5,4.5,0.004214929],[3,4.5,0.004759113],[3.5,4.5,0.005222199],[4,4.5,0.005568937],[4.5,4.5,0.005771422],[5,4.5,0.005812794],[5.5,4.5,0.005689559],[6,4.5,0.005412076],[6.5,4.5,0.005003118],[7,4.5,0.004494788],[7.5,4.5,0.003924364],[8,4.5,0.003329821],[8.5,4.5,0.002745769],[9,4.5,0.002200386],[9.5,4.5,0.001713662],[10,4.5,0.00129701],[10.5,4.5,0.0009540096],[11,4.5,0.0006819521],[11.5,4.5,0.0004737471],[12,4.5,0.0003198386],[12.5,4.5,0.000209849],[13,4.5,0.0001338056],[13.5,4.5,8.291505e-005],[14,4.5,4.993258e-005],[14.5,4.5,2.92231e-005],[15,4.5,1.662111e-005],[-5,5,1.200921e-005],[-4.5,5,2.096421e-005],[-4,5,3.556592e-005],[-3.5,5,5.86383e-005],[-3,5,9.395506e-005],[-2.5,5,0.0001463021],[-2,5,0.0002213974],[-1.5,5,0.0003256013],[-1,5,0.0004653625],[-0.5,5,0.0006463803],[0,5,0.0008725221],[0.5,5,0.001144607],[1,5,0.001459244],[1.5,5,0.001807969],[2,5,0.002176936],[2.5,5,0.00254737],[3,5,0.002896876],[3.5,5,0.003201543],[4,5,0.003438589],[4.5,5,0.003589161],[5,5,0.003640803],[5.5,5,0.003589161],[6,5,0.003438589],[6.5,5,0.003201543],[7,5,0.002896876],[7.5,5,0.00254737],[8,5,0.002176936],[8.5,5,0.001807969],[9,5,0.001459244],[9.5,5,0.001144607],[10,5,0.0008725221],[10.5,5,0.0006463803],[11,5,0.0004653625],[11.5,5,0.0003256013],[12,5,0.0002213974],[12.5,5,0.0001463021],[13,5,9.395506e-005],[13.5,5,5.86383e-005],[14,5,3.556592e-005],[14.5,5,2.096421e-005],[15,5,1.200921e-005],[-5,5.5,6.114566e-006],[-4.5,5.5,1.075058e-005],[-4,5.5,1.836917e-005],[-3.5,5.5,3.050274e-005],[-3,5.5,4.922433e-005],[-2.5,5.5,7.719912e-005],[-2,5.5,0.0001176621],[-1.5,5.5,0.0001742818],[-1,5.5,0.0002508762],[-0.5,5.5,0.0003509605],[0,5.5,0.0004771432],[0.5,5.5,0.0006304211],[1,5.5,0.0008094767],[1.5,5.5,0.001010112],[2,5.5,0.001224973],[2.5,5.5,0.001443693],[3,5.5,0.00165354],[3.5,5.5,0.001840544],[4,5.5,0.001990991],[4.5,5.5,0.002093072],[5,5.5,0.002138407],[5.5,5.5,0.002123188],[6,5.5,0.002048697],[6.5,5.5,0.001921139],[7,5.5,0.00175078],[7.5,5.5,0.001550586],[8,5.5,0.001334602],[8.5,5.5,0.001116347],[9,5.5,0.000907482],[9.5,5.5,0.0007169166],[10,5.5,0.0005504157],[10.5,5.5,0.000410681],[11,5.5,0.0002977898],[11.5,5.5,0.000209849],[12,5.5,0.0001437128],[12.5,5.5,9.56479e-005],[13,5.5,6.18653e-005],[13.5,5.5,3.888753e-005],[14,5.5,2.375556e-005],[14.5,5.5,1.410301e-005],[15,5.5,8.136727e-006],[-5,6,2.919429e-006],[-4.5,6,5.16971e-006],[-4,6,8.896641e-006],[-3.5,6,1.487913e-005],[-3,6,2.418358e-005],[-2.5,6,3.819928e-005],[-2,6,5.86383e-005],[-1.5,6,8.747806e-005],[-1,6,0.0001268261],[-0.5,6,0.0001786938],[0,6,0.000244682],[0.5,6,0.0003256013],[1,6,0.0004210774],[1.5,6,0.0005292114],[2,6,0.0006463803],[2.5,6,0.000767253],[3,6,0.0008850762],[3.5,6,0.0009922344],[4,6,0.001081035],[4.5,6,0.001144607],[5,6,0.001177782],[5.5,6,0.001177782],[6,6,0.001144607],[6.5,6,0.001081035],[7,6,0.0009922344],[7.5,6,0.0008850762],[8,6,0.000767253],[8.5,6,0.0006463803],[9,6,0.0005292114],[9.5,6,0.0004210774],[10,6,0.0003256013],[10.5,6,0.000244682],[11,6,0.0001786938],[11.5,6,0.0001268261],[12,6,8.747806e-005],[12.5,6,5.86383e-005],[13,6,3.819928e-005],[13.5,6,2.418358e-005],[14,6,1.487913e-005],[14.5,6,8.896641e-006],[15,6,5.16971e-006],[-5,6.5,1.307107e-006],[-4.5,6.5,2.331211e-006],[-4,6.5,4.040579e-006],[-3.5,6.5,6.806082e-006],[-3,6.5,1.114146e-005],[-2.5,6.5,1.772471e-005],[-2,6.5,2.740358e-005],[-1.5,6.5,4.11744e-005],[-1,6.5,6.012273e-005],[-0.5,6.5,8.531822e-005],[0,6.5,0.0001176621],[0.5,6.5,0.0001576967],[1,6.5,0.0002054],[1.5,6.5,0.000259998],[2,6.5,0.0003198386],[2.5,6.5,0.0003823697],[3,6.5,0.0004442502],[3.5,6.5,0.0005016068],[4,6.5,0.0005504157],[4.5,6.5,0.0005869616],[5,6.5,0.0006083034],[5.5,6.5,0.000612664],[6,6.5,0.0005996751],[6.5,6.5,0.0005704286],[7,6.5,0.0005273247],[7.5,6.5,0.0004737471],[8,6.5,0.0004136249],[8.5,6.5,0.0003509605],[9,6.5,0.000289402],[9.5,6.5,0.000231919],[10,6.5,0.0001806187],[10.5,6.5,0.0001367038],[11,6.5,0.0001005519],[11.5,6.5,7.187723e-005],[12,6.5,4.993258e-005],[12.5,6.5,3.371074e-005],[13,6.5,2.211792e-005],[13.5,6.5,1.410301e-005],[14,6.5,8.739182e-006],[14.5,6.5,5.262856e-006],[15,6.5,3.080092e-006],[-5,7,5.487894e-007],[-4.5,7,9.857758e-007],[-4,7,1.720847e-006],[-3.5,7,2.919429e-006],[-3,7,4.813325e-006],[-2.5,7,7.712301e-006],[-2,7,1.200921e-005],[-1.5,7,1.817341e-005],[-1,7,2.672699e-005],[-0.5,7,3.819928e-005],[0,7,5.305813e-005],[0.5,7,7.162098e-005],[1,7,9.395506e-005],[1.5,7,0.000119782],[2,7,0.0001484071],[2.5,7,0.0001786938],[3,7,0.0002091008],[3.5,7,0.00023779],[4,7,0.0002627986],[4.5,7,0.0002822566],[5,7,0.0002946163],[5.5,7,0.0002988553],[6,7,0.0002946163],[6.5,7,0.0002822566],[7,7,0.0002627986],[7.5,7,0.00023779],[8,7,0.0002091008],[8.5,7,0.0001786938],[9,7,0.0001484071],[9.5,7,0.000119782],[10,7,9.395506e-005],[10.5,7,7.162098e-005],[11,7,5.305813e-005],[11.5,7,3.819928e-005],[12,7,2.672699e-005],[12.5,7,1.817341e-005],[13,7,1.200921e-005],[13.5,7,7.712301e-006],[14,7,4.813325e-006],[14.5,7,2.919429e-006],[15,7,1.720847e-006],[-5,7.5,2.160634e-007],[-4.5,7.5,3.908911e-007],[-4,7.5,6.872614e-007],[-3.5,7.5,1.174302e-006],[-3,7.5,1.949975e-006],[-2.5,7.5,3.146806e-006],[-2,7.5,4.935174e-006],[-1.5,7.5,7.521884e-006],[-1,7.5,1.114146e-005],[-0.5,7.5,1.603798e-005],[0,7.5,2.243616e-005],[0.5,7.5,3.050274e-005],[1,7.5,4.030147e-005],[1.5,7.5,5.174812e-005],[2,7.5,6.457431e-005],[2.5,7.5,7.830988e-005],[3,7.5,9.229218e-005],[3.5,7.5,0.0001057073],[4,7.5,0.0001176621],[4.5,7.5,0.0001272798],[5,7.5,0.0001338056],[5.5,7.5,0.0001367038],[6,7.5,0.0001357308],[6.5,7.5,0.0001309688],[7,7.5,0.0001228143],[7.5,7.5,0.0001119236],[8,7.5,9.912563e-005],[8.5,7.5,8.531822e-005],[9,7.5,7.136565e-005],[9.5,7.5,5.801339e-005],[10,7.5,4.583094e-005],[10.5,7.5,3.518689e-005],[11,7.5,2.625395e-005],[11.5,7.5,1.903707e-005],[12,7.5,1.34152e-005],[12.5,7.5,9.18725e-006],[13,7.5,6.114566e-006],[13.5,7.5,3.954916e-006],[14,7.5,2.485997e-006],[14.5,7.5,1.518642e-006],[15,7.5,9.01575e-007],[-5,8,7.976966e-008],[-4.5,8,1.453498e-007],[-4,8,2.573847e-007],[-3.5,8,4.429376e-007],[-3,8,7.407882e-007],[-2.5,8,1.204029e-006],[-2,8,1.90183e-006],[-1.5,8,2.919429e-006],[-1,8,4.355276e-006],[-0.5,8,6.314298e-006],[0,8,8.896641e-006],[0.5,8,1.2182e-005],[1,8,1.621074e-005],[1.5,8,2.096421e-005],[2,8,2.634789e-005],[2.5,8,3.218138e-005],[3,8,3.819928e-005],[3.5,8,4.406535e-005],[4,8,4.940045e-005],[4.5,8,5.382154e-005],[5,8,5.698662e-005],[5.5,8,5.86383e-005],[6,8,5.86383e-005],[6.5,8,5.698662e-005],[7,8,5.382154e-005],[7.5,8,4.940045e-005],[8,8,4.406535e-005],[8.5,8,3.819928e-005],[9,8,3.218138e-005],[9.5,8,2.634789e-005],[10,8,2.096421e-005],[10.5,8,1.621074e-005],[11,8,1.2182e-005],[11.5,8,8.896641e-006],[12,8,6.314298e-006],[12.5,8,4.355276e-006],[13,8,2.919429e-006],[13.5,8,1.90183e-006],[14,8,1.204029e-006],[14.5,8,7.407882e-007],[15,8,4.429376e-007],[-5,8.5,2.761693e-008],[-4.5,8.5,5.068205e-008],[-4,8.5,9.039088e-008],[-3.5,8.5,1.566703e-007],[-3,8.5,2.639004e-007],[-2.5,8.5,4.320015e-007],[-2,8.5,6.872614e-007],[-1.5,8.5,1.062552e-006],[-1,8.5,1.596505e-006],[-0.5,8.5,2.331211e-006],[0,8.5,3.308146e-006],[0.5,8.5,4.562253e-006],[1,8.5,6.114566e-006],[1.5,8.5,7.964223e-006],[2,8.5,1.008121e-005],[2.5,8.5,1.240149e-005],[3,8.5,1.482608e-005],[3.5,8.5,1.722545e-005],[4,8.5,1.944941e-005],[4.5,8.5,2.134193e-005],[5,8.5,2.275897e-005],[5.5,8.5,2.358648e-005],[6,8.5,2.375556e-005],[6.5,8.5,2.325193e-005],[7,8.5,2.211792e-005],[7.5,8.5,2.04466e-005],[8,8.5,1.836917e-005],[8.5,8.5,1.603798e-005],[9,8.5,1.360822e-005],[9.5,8.5,1.122133e-005],[10,8.5,8.992474e-006],[10.5,8.5,7.003346e-006],[11,8.5,5.300582e-006],[11.5,8.5,3.898819e-006],[12,8.5,2.786983e-006],[12.5,8.5,1.936096e-006],[13,8.5,1.307107e-006],[13.5,8.5,8.576047e-007],[14,8.5,5.468329e-007],[14.5,8.5,3.388549e-007],[15,8.5,2.040631e-007],[-5,9,8.965904e-009],[-4.5,9,1.6572e-008],[-4,9,2.976785e-008],[-3.5,9,5.196508e-008],[-3,9,8.815911e-008],[-2.5,9,1.453498e-007],[-2,9,2.328913e-007],[-1.5,9,3.626467e-007],[-1,9,5.487894e-007],[-0.5,9,8.07085e-007],[0,9,1.153518e-006],[0.5,9,1.602217e-006],[1,9,2.162766e-006],[1.5,9,2.837197e-006],[2,9,3.617104e-006],[2.5,9,4.481507e-006],[3,9,5.396085e-006],[3.5,9,6.314298e-006],[4,9,7.180637e-006],[4.5,9,7.935831e-006],[5,9,8.523411e-006],[5.5,9,8.896641e-006],[6,9,9.024648e-006],[6.5,9,8.896641e-006],[7,9,8.523411e-006],[7.5,9,7.935831e-006],[8,9,7.180637e-006],[8.5,9,6.314298e-006],[9,9,5.396085e-006],[9.5,9,4.481507e-006],[10,9,3.617104e-006],[10.5,9,2.837197e-006],[11,9,2.162766e-006],[11.5,9,1.602217e-006],[12,9,1.153518e-006],[12.5,9,8.07085e-007],[13,9,5.487894e-007],[13.5,9,3.626467e-007],[14,9,2.328913e-007],[14.5,9,1.453498e-007],[15,9,8.815911e-008],[-5,9.5,2.729568e-009],[-4.5,9.5,5.081324e-009],[-4,9.5,9.192878e-009],[-3.5,9.5,1.616284e-008],[-3,9.5,2.761693e-008],[-2.5,9.5,4.585902e-008],[-2,9.5,7.400579e-008],[-1.5,9.5,1.160642e-007],[-1,9.5,1.768978e-007],[-0.5,9.5,2.620222e-007],[0,9.5,3.771771e-007],[0.5,9.5,5.276478e-007],[1,9.5,7.173558e-007],[1.5,9.5,9.477998e-007],[2,9.5,1.216999e-006],[2.5,9.5,1.518642e-006],[3,9.5,1.841672e-006],[3.5,9.5,2.170504e-006],[4,9.5,2.485997e-006],[4.5,9.5,2.767146e-006],[5,9.5,2.993335e-006],[5.5,9.5,3.146806e-006],[6,9.5,3.214965e-006],[6.5,9.5,3.192083e-006],[7,9.5,3.080092e-006],[7.5,9.5,2.888316e-006],[8,9.5,2.632191e-006],[8.5,9.5,2.331211e-006],[9,9.5,2.006492e-006],[9.5,9.5,1.678359e-006],[10,9.5,1.364344e-006],[10.5,9.5,1.07784e-006],[11,9.5,8.275165e-007],[11.5,9.5,6.174338e-007],[12,9.5,4.477089e-007],[12.5,9.5,3.154952e-007],[13,9.5,2.160634e-007],[13.5,9.5,1.438008e-007],[14,9.5,9.301072e-008],[14.5,9.5,5.846504e-008],[15,9.5,3.571504e-008],[-5,10,7.792463e-010],[-4.5,10,1.461032e-009],[-4,10,2.662175e-009],[-3.5,10,4.714166e-009],[-3,10,8.112687e-009],[-2.5,10,1.356801e-008],[-2,10,2.205257e-008],[-1.5,10,3.483323e-008],[-1,10,5.347121e-008],[-0.5,10,7.976966e-008],[0,10,1.156504e-007],[0.5,10,1.629477e-007],[1,10,2.231211e-007],[1.5,10,2.9691e-007],[2,10,3.839729e-007],[2.5,10,4.825783e-007],[3,10,5.894225e-007],[3.5,10,6.996441e-007],[4,10,8.07085e-007],[4.5,10,9.048007e-007],[5,10,9.857758e-007],[5.5,10,1.043746e-006],[6,10,1.073998e-006],[6.5,10,1.073998e-006],[7,10,1.043746e-006],[7.5,10,9.857758e-007],[8,10,9.048007e-007],[8.5,10,8.07085e-007],[9,10,6.996441e-007],[9.5,10,5.894225e-007],[10,10,4.825783e-007],[10.5,10,3.839729e-007],[11,10,2.9691e-007],[11.5,10,2.231211e-007],[12,10,1.629477e-007],[12.5,10,1.156504e-007],[13,10,7.976966e-008],[13.5,10,5.347121e-008],[14,10,3.483323e-008],[14.5,10,2.205257e-008],[15,10,1.356801e-008],[-5,10.5,2.086107e-010],[-4.5,10.5,3.939344e-010],[-4,10.5,7.229407e-010],[-3.5,10.5,1.289357e-009],[-3,10.5,2.234781e-009],[-2.5,10.5,3.764337e-009],[-2,10.5,6.16217e-009],[-1.5,10.5,9.803258e-009],[-1,10.5,1.515649e-008],[-0.5,10.5,2.277291e-008],[0,10.5,3.325295e-008],[0.5,10.5,4.718817e-008],[1,10.5,6.507705e-008],[1.5,10.5,8.72196e-008],[2,10.5,1.136035e-007],[2.5,10.5,1.438008e-007],[3,10.5,1.768978e-007],[3.5,10.5,2.114827e-007],[4,10.5,2.457079e-007],[4.5,10.5,2.774309e-007],[5,10.5,3.044263e-007],[5.5,10.5,3.246393e-007],[6,10.5,3.364431e-007],[6.5,10.5,3.388549e-007],[7,10.5,3.316709e-007],[7.5,10.5,3.154952e-007],[8,10.5,2.916551e-007],[8.5,10.5,2.620222e-007],[9,10.5,2.287695e-007],[9.5,10.5,1.941108e-007],[10,10.5,1.600637e-007],[10.5,10.5,1.282708e-007],[11,10.5,9.989736e-008],[11.5,10.5,7.560874e-008],[12,10.5,5.561367e-008],[12.5,10.5,3.975417e-008],[13,10.5,2.761693e-008],[13.5,10.5,1.864489e-008],[14,10.5,1.223307e-008],[14.5,10.5,7.800153e-009],[15,10.5,4.833505e-009],[-5,11,5.236965e-011],[-4.5,11,9.960221e-011],[-4,11,1.840983e-010],[-3.5,11,3.30691e-010],[-3,11,5.772798e-010],[-2.5,11,9.793593e-010],[-2,11,1.61469e-009],[-1.5,11,2.587189e-009],[-1,11,4.028641e-009],[-0.5,11,6.096499e-009],[0,11,8.965904e-009],[0.5,11,1.281443e-008],[1,11,1.779902e-008],[1.5,11,2.402616e-008],[2,11,3.151841e-008],[2.5,11,4.01824e-008],[3,11,4.978505e-008],[3.5,11,5.994509e-008],[4,11,7.014552e-008],[4.5,11,7.976966e-008],[5,11,8.815911e-008],[5.5,11,9.468654e-008],[6,11,9.883275e-008],[6.5,11,1.002548e-007],[7,11,9.883275e-008],[7.5,11,9.468654e-008],[8,11,8.815911e-008],[8.5,11,7.976966e-008],[9,11,7.014552e-008],[9.5,11,5.994509e-008],[10,11,4.978505e-008],[10.5,11,4.01824e-008],[11,11,3.151841e-008],[11.5,11,2.402616e-008],[12,11,1.779902e-008],[12.5,11,1.281443e-008],[13,11,8.965904e-009],[13.5,11,6.096499e-009],[14,11,4.028641e-009],[14.5,11,2.587189e-009],[15,11,1.61469e-009],[-5,11.5,1.232831e-011],[-4.5,11.5,2.361539e-011],[-4,11.5,4.396207e-011],[-3.5,11.5,7.953398e-011],[-3,11.5,1.39836e-010],[-2.5,11.5,2.389333e-010],[-2,11.5,3.967583e-010],[-1.5,11.5,6.402755e-010],[-1,11.5,1.004152e-009],[-0.5,11.5,1.530465e-009],[0,11.5,2.266936e-009],[0.5,11.5,3.263221e-009],[1,11.5,4.565048e-009],[1.5,11.5,6.206343e-009],[2,11.5,8.200075e-009],[2.5,11.5,1.05291e-008],[3,11.5,1.313883e-008],[3.5,11.5,1.593358e-008],[4,11.5,1.877854e-008],[4.5,11.5,2.150809e-008],[5,11.5,2.394051e-008],[5.5,11.5,2.589742e-008],[6,11.5,2.722521e-008],[6.5,11.5,2.78149e-008],[7,11.5,2.761693e-008],[7.5,11.5,2.664802e-008],[8,11.5,2.498883e-008],[8.5,11.5,2.277291e-008],[9,11.5,2.016893e-008],[9.5,11.5,1.735956e-008],[10,11.5,1.452065e-008],[10.5,11.5,1.180389e-008],[11,11.5,9.325147e-009],[11.5,11.5,7.15942e-009],[12,11.5,5.341849e-009],[12.5,11.5,3.873441e-009],[13,11.5,2.729568e-009],[13.5,11.5,1.869315e-009],[14,11.5,1.244121e-009],[14.5,11.5,8.047004e-010],[15,11.5,5.058217e-010],[-5,12,2.721503e-012],[-4.5,12,5.25052e-012],[-4,12,9.844359e-012],[-3.5,12,1.793759e-011],[-3,12,3.17638e-011],[-2.5,12,5.466285e-011],[-2,12,9.142052e-011],[-1.5,12,1.48589e-010],[-1,12,2.347045e-010],[-0.5,12,3.602862e-010],[0,12,5.374838e-010],[0.5,12,7.792463e-010],[1,12,1.097933e-009],[1.5,12,1.503378e-009],[2,12,2.000564e-009],[2.5,12,2.587189e-009],[3,12,3.251587e-009],[3.5,12,3.971498e-009],[4,12,4.714166e-009],[4.5,12,5.438096e-009],[5,12,6.096499e-009],[5.5,12,6.642106e-009],[6,12,7.032708e-009],[6.5,12,7.236541e-009],[7,12,7.236541e-009],[7.5,12,7.032708e-009],[8,12,6.642106e-009],[8.5,12,6.096499e-009],[9,12,5.438096e-009],[9.5,12,4.714166e-009],[10,12,3.971498e-009],[10.5,12,3.251587e-009],[11,12,2.587189e-009],[11.5,12,2.000564e-009],[12,12,1.503378e-009],[12.5,12,1.097933e-009],[13,12,7.792463e-010],[13.5,12,5.374838e-010],[14,12,3.602862e-010],[14.5,12,2.347045e-010],[15,12,1.48589e-010],[-5,12.5,5.633717e-013],[-4.5,12.5,1.094689e-012],[-4,12.5,2.067178e-012],[-3.5,12.5,3.793645e-012],[-3,12.5,6.765922e-012],[-2.5,12.5,1.172706e-011],[-2,12.5,1.975343e-011],[-1.5,12.5,3.23361e-011],[-1,12.5,5.144277e-011],[-0.5,12.5,7.953398e-011],[0,12.5,1.195013e-010],[0.5,12.5,1.744954e-010],[1,12.5,2.476208e-010],[1.5,12.5,3.41493e-010],[2,12.5,4.576864e-010],[2.5,12.5,5.961367e-010],[3,12.5,7.545973e-010],[3.5,12.5,9.282742e-010],[4,12.5,1.10976e-009],[4.5,12.5,1.289357e-009],[5,12.5,1.455824e-009],[5.5,12.5,1.597483e-009],[6,12.5,1.703551e-009],[6.5,12.5,1.765491e-009],[7,12.5,1.778147e-009],[7.5,12.5,1.740449e-009],[8,12.5,1.655567e-009],[8.5,12.5,1.530465e-009],[9,12.5,1.374966e-009],[9.5,12.5,1.200472e-009],[10,12.5,1.0186e-009],[10.5,12.5,8.399372e-010],[11,12.5,6.731031e-010],[11.5,12.5,5.242133e-010],[12,12.5,3.967583e-010],[12.5,12.5,2.918338e-010],[13,12.5,2.086107e-010],[13.5,12.5,1.449204e-010],[14,12.5,9.783938e-011],[14.5,12.5,6.419328e-011],[15,12.5,4.093145e-011],[-5,13,1.093609e-013],[-4.5,13,2.140228e-013],[-4,13,4.070515e-013],[-3.5,13,7.523678e-013],[-3,13,1.351458e-012],[-2.5,13,2.359211e-012],[-2,13,4.002418e-012],[-1.5,13,6.598871e-012],[-1,13,1.057325e-011],[-0.5,13,1.646413e-011],[0,13,2.4915e-011],[0.5,13,3.66416e-011],[1,13,5.236965e-011],[1.5,13,7.274052e-011],[2,13,9.818943e-011],[2.5,13,1.288085e-010],[3,13,1.642163e-010],[3.5,13,2.034601e-010],[4,13,2.449819e-010],[4.5,13,2.866687e-010],[5,13,3.260004e-010],[5.5,13,3.602862e-010],[6,13,3.869623e-010],[6.5,13,4.039069e-010],[7,13,4.097184e-010],[7.5,13,4.039069e-010],[8,13,3.869623e-010],[8.5,13,3.602862e-010],[9,13,3.260004e-010],[9.5,13,2.866687e-010],[10,13,2.449819e-010],[10.5,13,2.034601e-010],[11,13,1.642163e-010],[11.5,13,1.288085e-010],[12,13,9.818943e-011],[12.5,13,7.274052e-011],[13,13,5.236965e-011],[13.5,13,3.66416e-011],[14,13,2.4915e-011],[14.5,13,1.646413e-011],[15,13,1.057325e-011],[-5,13.5,1.990722e-014],[-4.5,13.5,3.923832e-014],[-4,13.5,7.516261e-014],[-3.5,13.5,1.399216e-013],[-3,13.5,2.531392e-013],[-2.5,13.5,4.450672e-013],[-2,13.5,7.604722e-013],[-1.5,13.5,1.262795e-012],[-1,13.5,2.037857e-012],[-0.5,13.5,3.195995e-012],[0,13.5,4.871136e-012],[0.5,13.5,7.215159e-012],[1,13.5,1.038612e-011],[1.5,13.5,1.452955e-011],[2,13.5,1.975343e-011],[2.5,13.5,2.609904e-011],[3,13.5,3.351183e-011],[3.5,13.5,4.181801e-011],[4,13.5,5.07131e-011],[4.5,13.5,5.976798e-011],[5,13.5,6.845553e-011],[5.5,13.5,7.619739e-011],[6,13.5,8.242581e-011],[6.5,13.5,8.665187e-011],[7,13.5,8.852873e-011],[7.5,13.5,8.789864e-011],[8,13.5,8.481479e-011],[8.5,13.5,7.953398e-011],[9,13.5,7.248119e-011],[9.5,13.5,6.419328e-011],[10,13.5,5.525167e-011],[10.5,13.5,4.621605e-011],[11,13.5,3.756919e-011],[11.5,13.5,2.96799e-011],[12,13.5,2.278687e-011],[12.5,13.5,1.700193e-011],[13,13.5,1.232831e-011],[13.5,13.5,8.687617e-012],[14,13.5,5.949619e-012],[14.5,13.5,3.959764e-012],[15,13.5,2.561185e-012],[-5,14,3.39813e-015],[-4.5,14,6.745931e-015],[-4,14,1.301474e-014],[-3.5,14,2.440173e-014],[-3,14,4.446284e-014],[-2.5,14,7.873458e-014],[-2,14,1.354957e-013],[-1.5,14,2.266088e-013],[-1,14,3.683154e-013],[-0.5,14,5.817742e-013],[0,14,8.930601e-013],[0.5,14,1.332289e-012],[1,14,1.931558e-012],[1.5,14,2.721503e-012],[2,14,3.726502e-012],[2.5,14,4.958902e-012],[3,14,6.413e-012],[3.5,14,8.059879e-012],[4,14,9.844359e-012],[4.5,14,1.168525e-011],[5,14,1.347969e-011],[5.5,14,1.511171e-011],[6,14,1.646413e-011],[6.5,14,1.743234e-011],[7,14,1.793759e-011],[7.5,14,1.793759e-011],[8,14,1.743234e-011],[8.5,14,1.646413e-011],[9,14,1.511171e-011],[9.5,14,1.347969e-011],[10,14,1.168525e-011],[10.5,14,9.844359e-012],[11,14,8.059879e-012],[11.5,14,6.413e-012],[12,14,4.958902e-012],[12.5,14,3.726502e-012],[13,14,2.721503e-012],[13.5,14,1.931558e-012],[14,14,1.332289e-012],[14.5,14,8.930601e-013],[15,14,5.817742e-013],[-5,14.5,5.439393e-016],[-4.5,14.5,1.087563e-015],[-4,14.5,2.113246e-015],[-3.5,14.5,3.990592e-015],[-3,14.5,7.323457e-015],[-2.5,14.5,1.30613e-014],[-2,14.5,2.263854e-014],[-1.5,14.5,3.813309e-014],[-1,14.5,6.242336e-014],[-0.5,14.5,9.930792e-014],[0,14.5,1.535367e-013],[0.5,14.5,2.306918e-013],[1,14.5,3.368554e-013],[1.5,14.5,4.780206e-013],[2,14.5,6.592366e-013],[2.5,14.5,8.835426e-013],[3,14.5,1.150815e-012],[3.5,14.5,1.456716e-012],[4,14.5,1.791991e-012],[4.5,14.5,2.14234e-012],[5,14.5,2.489044e-012],[5.5,14.5,2.810401e-012],[6,14.5,3.083867e-012],[6.5,14.5,3.288627e-012],[7,14.5,3.4082e-012],[7.5,14.5,3.432632e-012],[8,14.5,3.359858e-012],[8.5,14.5,3.195995e-012],[9,14.5,2.954493e-012],[9.5,14.5,2.654309e-012],[10,14.5,2.317456e-012],[10.5,14.5,1.96636e-012],[11,14.5,1.62146e-012],[11.5,14.5,1.299395e-012],[12,14.5,1.01197e-012],[12.5,14.5,7.659236e-013],[13,14.5,5.633717e-013],[13.5,14.5,4.027135e-013],[14,14.5,2.797621e-013],[14.5,14.5,1.888744e-013],[15,14.5,1.239222e-013],[-5,15,8.164733e-017],[-4.5,15,1.644175e-016],[-4,15,3.217702e-016],[-3.5,15,6.119771e-016],[-3,15,1.131139e-015],[-2.5,15,2.031836e-015],[-2,15,3.54693e-015],[-1.5,15,6.017392e-015],[-1,15,9.921001e-015],[-0.5,15,1.589624e-014],[0,15,2.475282e-014],[0.5,15,3.745819e-014],[1,15,5.508842e-014],[1.5,15,7.873458e-014],[2,15,1.093609e-013],[2.5,15,1.476218e-013],[3,15,1.936558e-013],[3.5,15,2.468892e-013],[4,15,3.058899e-013],[4.5,15,3.683154e-013],[5,15,4.30989e-013],[5.5,15,4.901218e-013],[6,15,5.416683e-013],[6.5,15,5.817742e-013],[7,15,6.072494e-013],[7.5,15,6.159866e-013],[8,15,6.072494e-013],[8.5,15,5.817742e-013],[9,15,5.416683e-013],[9.5,15,4.901218e-013],[10,15,4.30989e-013],[10.5,15,3.683154e-013],[11,15,3.058899e-013],[11.5,15,2.468892e-013],[12,15,1.936558e-013],[12.5,15,1.476218e-013],[13,15,1.093609e-013],[13.5,15,7.873458e-014],[14,15,5.508842e-014],[14.5,15,3.745819e-014],[15,15,2.475282e-014]],"colors":[[0,1,0,1]],"centers":[[-5,-5,1.200921e-005],[-4.5,-5,1.817341e-005],[-4,-5,2.672699e-005],[-3.5,-5,3.819928e-005],[-3,-5,5.305813e-005],[-2.5,-5,7.162098e-005],[-2,-5,9.395506e-005],[-1.5,-5,0.000119782],[-1,-5,0.0001484071],[-0.5,-5,0.0001786938],[0,-5,0.0002091008],[0.5,-5,0.00023779],[1,-5,0.0002627986],[1.5,-5,0.0002822566],[2,-5,0.0002946163],[2.5,-5,0.0002988553],[3,-5,0.0002946163],[3.5,-5,0.0002822566],[4,-5,0.0002627986],[4.5,-5,0.00023779],[5,-5,0.0002091008],[5.5,-5,0.0001786938],[6,-5,0.0001484071],[6.5,-5,0.000119782],[7,-5,9.395506e-005],[7.5,-5,7.162098e-005],[8,-5,5.305813e-005],[8.5,-5,3.819928e-005],[9,-5,2.672699e-005],[9.5,-5,1.817341e-005],[10,-5,1.200921e-005],[10.5,-5,7.712301e-006],[11,-5,4.813325e-006],[11.5,-5,2.919429e-006],[12,-5,1.720847e-006],[12.5,-5,9.857758e-007],[13,-5,5.487894e-007],[13.5,-5,2.9691e-007],[14,-5,1.561117e-007],[14.5,-5,7.976966e-008],[15,-5,3.961244e-008],[-5,-4.5,2.211792e-005],[-4.5,-4.5,3.371074e-005],[-4,-4.5,4.993258e-005],[-3.5,-4.5,7.187723e-005],[-3,-4.5,0.0001005519],[-2.5,-4.5,0.0001367038],[-2,-4.5,0.0001806187],[-1.5,-4.5,0.000231919],[-1,-4.5,0.000289402],[-0.5,-4.5,0.0003509605],[0,-4.5,0.0004136249],[0.5,-4.5,0.0004737471],[1,-4.5,0.0005273247],[1.5,-4.5,0.0005704286],[2,-4.5,0.0005996751],[2.5,-4.5,0.000612664],[3,-4.5,0.0006083034],[3.5,-4.5,0.0005869616],[4,-4.5,0.0005504157],[4.5,-4.5,0.0005016068],[5,-4.5,0.0004442502],[5.5,-4.5,0.0003823697],[6,-4.5,0.0003198386],[6.5,-4.5,0.000259998],[7,-4.5,0.0002054],[7.5,-4.5,0.0001576967],[8,-4.5,0.0001176621],[8.5,-4.5,8.531822e-005],[9,-4.5,6.012273e-005],[9.5,-4.5,4.11744e-005],[10,-4.5,2.740358e-005],[10.5,-4.5,1.772471e-005],[11,-4.5,1.114146e-005],[11.5,-4.5,6.806082e-006],[12,-4.5,4.040579e-006],[12.5,-4.5,2.331211e-006],[13,-4.5,1.307107e-006],[13.5,-4.5,7.1225e-007],[14,-4.5,3.771771e-007],[14.5,-4.5,1.941108e-007],[15,-4.5,9.708354e-008],[-5,-4,3.819928e-005],[-4.5,-4,5.86383e-005],[-4,-4,8.747806e-005],[-3.5,-4,0.0001268261],[-3,-4,0.0001786938],[-2.5,-4,0.000244682],[-2,-4,0.0003256013],[-1.5,-4,0.0004210774],[-1,-4,0.0005292114],[-0.5,-4,0.0006463803],[0,-4,0.000767253],[0.5,-4,0.0008850762],[1,-4,0.0009922344],[1.5,-4,0.001081035],[2,-4,0.001144607],[2.5,-4,0.001177782],[3,-4,0.001177782],[3.5,-4,0.001144607],[4,-4,0.001081035],[4.5,-4,0.0009922344],[5,-4,0.0008850762],[5.5,-4,0.000767253],[6,-4,0.0006463803],[6.5,-4,0.0005292114],[7,-4,0.0004210774],[7.5,-4,0.0003256013],[8,-4,0.000244682],[8.5,-4,0.0001786938],[9,-4,0.0001268261],[9.5,-4,8.747806e-005],[10,-4,5.86383e-005],[10.5,-4,3.819928e-005],[11,-4,2.418358e-005],[11.5,-4,1.487913e-005],[12,-4,8.896641e-006],[12.5,-4,5.16971e-006],[13,-4,2.919429e-006],[13.5,-4,1.602217e-006],[14,-4,8.545474e-007],[14.5,-4,4.429376e-007],[15,-4,2.231211e-007],[-5,-3.5,6.18653e-005],[-4.5,-3.5,9.56479e-005],[-4,-3.5,0.0001437128],[-3.5,-3.5,0.000209849],[-3,-3.5,0.0002977898],[-2.5,-3.5,0.000410681],[-2,-3.5,0.0005504157],[-1.5,-3.5,0.0007169166],[-1,-3.5,0.000907482],[-0.5,-3.5,0.001116347],[0,-3.5,0.001334602],[0.5,-3.5,0.001550586],[1,-3.5,0.00175078],[1.5,-3.5,0.001921139],[2,-3.5,0.002048697],[2.5,-3.5,0.002123188],[3,-3.5,0.002138407],[3.5,-3.5,0.002093072],[4,-3.5,0.001990991],[4.5,-3.5,0.001840544],[5,-3.5,0.00165354],[5.5,-3.5,0.001443693],[6,-3.5,0.001224973],[6.5,-3.5,0.001010112],[7,-3.5,0.0008094767],[7.5,-3.5,0.0006304211],[8,-3.5,0.0004771432],[8.5,-3.5,0.0003509605],[9,-3.5,0.0002508762],[9.5,-3.5,0.0001742818],[10,-3.5,0.0001176621],[10.5,-3.5,7.719912e-005],[11,-3.5,4.922433e-005],[11.5,-3.5,3.050274e-005],[12,-3.5,1.836917e-005],[12.5,-3.5,1.075058e-005],[13,-3.5,6.114566e-006],[13.5,-3.5,3.3798e-006],[14,-3.5,1.815549e-006],[14.5,-3.5,9.477998e-007],[15,-3.5,4.808579e-007],[-5,-3,9.395506e-005],[-4.5,-3,0.0001463021],[-4,-3,0.0002213974],[-3.5,-3,0.0003256013],[-3,-3,0.0004653625],[-2.5,-3,0.0006463803],[-2,-3,0.0008725221],[-1.5,-3,0.001144607],[-1,-3,0.001459244],[-0.5,-3,0.001807969],[0,-3,0.002176936],[0.5,-3,0.00254737],[1,-3,0.002896876],[1.5,-3,0.003201543],[2,-3,0.003438589],[2.5,-3,0.003589161],[3,-3,0.003640803],[3.5,-3,0.003589161],[4,-3,0.003438589],[4.5,-3,0.003201543],[5,-3,0.002896876],[5.5,-3,0.00254737],[6,-3,0.002176936],[6.5,-3,0.001807969],[7,-3,0.001459244],[7.5,-3,0.001144607],[8,-3,0.0008725221],[8.5,-3,0.0006463803],[9,-3,0.0004653625],[9.5,-3,0.0003256013],[10,-3,0.0002213974],[10.5,-3,0.0001463021],[11,-3,9.395506e-005],[11.5,-3,5.86383e-005],[12,-3,3.556592e-005],[12.5,-3,2.096421e-005],[13,-3,1.200921e-005],[13.5,-3,6.685623e-006],[14,-3,3.617104e-006],[14.5,-3,1.90183e-006],[15,-3,9.717935e-007],[-5,-2.5,0.0001338056],[-4.5,-2.5,0.000209849],[-4,-2.5,0.0003198386],[-3.5,-2.5,0.0004737471],[-3,-2.5,0.0006819521],[-2.5,-2.5,0.0009540096],[-2,-2.5,0.00129701],[-1.5,-2.5,0.001713662],[-1,-2.5,0.002200386],[-0.5,-2.5,0.002745769],[0,-2.5,0.003329821],[0.5,-2.5,0.003924364],[1,-2.5,0.004494788],[1.5,-2.5,0.005003118],[2,-2.5,0.005412076],[2.5,-2.5,0.005689559],[3,-2.5,0.005812794],[3.5,-2.5,0.005771422],[4,-2.5,0.005568937],[4.5,-2.5,0.005222199],[5,-2.5,0.004759113],[5.5,-2.5,0.004214929],[6,-2.5,0.003627823],[6.5,-2.5,0.003034544],[7,-2.5,0.002466792],[7.5,-2.5,0.001948781],[8,-2.5,0.001496185],[8.5,-2.5,0.001116347],[9,-2.5,0.0008094767],[9.5,-2.5,0.0005704286],[10,-2.5,0.0003906518],[10.5,-2.5,0.000259998],[11,-2.5,0.0001681673],[11.5,-2.5,0.0001057073],[12,-2.5,6.457431e-005],[12.5,-2.5,3.833595e-005],[13,-2.5,2.211792e-005],[13.5,-2.5,1.240149e-005],[14,-2.5,6.75764e-006],[14.5,-2.5,3.578556e-006],[15,-2.5,1.841672e-006],[-5,-2,0.0001786938],[-4.5,-2,0.0002822566],[-4,-2,0.0004332817],[-3.5,-2,0.0006463803],[-3,-2,0.0009371249],[-2.5,-2,0.001320379],[-2,-2,0.001807969],[-1.5,-2,0.002405887],[-1,-2,0.003111364],[-0.5,-2,0.003910373],[0,-2,0.00477614],[0.5,-2,0.005669275],[1,-2,0.006539878],[1.5,-2,0.007331676],[2,-2,0.007987824],[2.5,-2,0.008457565],[3,-2,0.008702695],[3.5,-2,0.008702695],[4,-2,0.008457565],[4.5,-2,0.007987824],[5,-2,0.007331676],[5.5,-2,0.006539878],[6,-2,0.005669275],[6.5,-2,0.00477614],[7,-2,0.003910373],[7.5,-2,0.003111364],[8,-2,0.002405887],[8.5,-2,0.001807969],[9,-2,0.001320379],[9.5,-2,0.0009371249],[10,-2,0.0006463803],[10.5,-2,0.0004332817],[11,-2,0.0002822566],[11.5,-2,0.0001786938],[12,-2,0.0001099427],[12.5,-2,6.573778e-005],[13,-2,3.819928e-005],[13.5,-2,2.157182e-005],[14,-2,1.183887e-005],[14.5,-2,6.314298e-006],[15,-2,3.272891e-006],[-5,-1.5,0.0002237823],[-4.5,-1.5,0.0003560102],[-4,-1.5,0.0005504157],[-3.5,-1.5,0.0008270098],[-3,-1.5,0.001207597],[-2.5,-1.5,0.001713662],[-2,-1.5,0.002363306],[-1.5,-1.5,0.003167423],[-1,-1.5,0.00412557],[-0.5,-1.5,0.005222199],[0,-1.5,0.006424131],[0.5,-1.5,0.007680101],[1,-1.5,0.008923004],[1.5,-1.5,0.01007504],[2,-1.5,0.01105539],[2.5,-1.5,0.01178944],[3,-1.5,0.0122181],[3.5,-1.5,0.01230568],[4,-1.5,0.0120448],[4.5,-1.5,0.01145737],[5,-1.5,0.0105916],[5.5,-1.5,0.009515466],[6,-1.5,0.008307878],[6.5,-1.5,0.00704923],[7,-1.5,0.005812794],[7.5,-1.5,0.004658217],[8,-1.5,0.003627823],[8.5,-1.5,0.002745769],[9,-1.5,0.002019638],[9.5,-1.5,0.001443693],[10,-1.5,0.001002923],[10.5,-1.5,0.0006770984],[11,-1.5,0.0004442502],[11.5,-1.5,0.0002832665],[12,-1.5,0.0001755312],[12.5,-1.5,0.0001057073],[13,-1.5,6.18653e-005],[13.5,-1.5,3.518689e-005],[14,-1.5,1.944941e-005],[14.5,-1.5,1.044776e-005],[15,-1.5,5.454212e-006],[-5,-1,0.0002627986],[-4.5,-1,0.0004210774],[-4,-1,0.0006556806],[-3.5,-1,0.0009922344],[-3,-1,0.001459244],[-2.5,-1,0.00208561],[-2,-1,0.002896876],[-1.5,-1,0.003910373],[-1,-1,0.005129773],[-0.5,-1,0.006539878],[0,-1,0.008102755],[0.5,-1,0.009756351],[1,-1,0.01141652],[1.5,-1,0.0129829],[2,-1,0.01434832],[2.5,-1,0.01541069],[3,-1,0.0160855],[3.5,-1,0.01631695],[4,-1,0.0160855],[4.5,-1,0.01541069],[5,-1,0.01434832],[5.5,-1,0.0129829],[6,-1,0.01141652],[6.5,-1,0.009756351],[7,-1,0.008102755],[7.5,-1,0.006539878],[8,-1,0.005129773],[8.5,-1,0.003910373],[9,-1,0.002896876],[9.5,-1,0.00208561],[10,-1,0.001459244],[10.5,-1,0.0009922344],[11,-1,0.0006556806],[11.5,-1,0.0004210774],[12,-1,0.0002627986],[12.5,-1,0.0001593954],[13,-1,9.395506e-005],[13.5,-1,5.382154e-005],[14,-1,2.996288e-005],[14.5,-1,1.621074e-005],[15,-1,8.523411e-006],[-5,-0.5,0.000289402],[-4.5,-0.5,0.0004670275],[-4,-0.5,0.0007324449],[-3.5,-0.5,0.001116347],[-3,-0.5,0.00165354],[-2.5,-0.5,0.002380247],[-2,-0.5,0.003329821],[-1.5,-0.5,0.004527009],[-1,-0.5,0.005981269],[-0.5,-0.5,0.007680101],[0,-0.5,0.009583677],[0.5,-0.5,0.01162222],[1,-0.5,0.01369737],[1.5,-0.5,0.01568835],[2,-0.5,0.0174626],[2.5,-0.5,0.01889],[3,-0.5,0.01985851],[3.5,-0.5,0.02028864],[4,-0.5,0.02014424],[4.5,-0.5,0.0194375],[5,-0.5,0.01822726],[5.5,-0.5,0.01661094],[6,-0.5,0.01471155],[6.5,-0.5,0.01266235],[7,-0.5,0.0105916],[7.5,-0.5,0.00860995],[8,-0.5,0.006801915],[8.5,-0.5,0.005222199],[9,-0.5,0.003896432],[9.5,-0.5,0.002825351],[10,-0.5,0.001990991],[10.5,-0.5,0.001363509],[11,-0.5,0.000907482],[11.5,-0.5,0.0005869616],[12,-0.5,0.0003689546],[12.5,-0.5,0.0002253865],[13,-0.5,0.0001338056],[13.5,-0.5,7.719912e-005],[14,-0.5,4.328545e-005],[14.5,-0.5,2.358648e-005],[15,-0.5,1.249039e-005],[-5,0,0.0002988553],[-4.5,0,0.0004857401],[-4,0,0.000767253],[-3.5,0,0.001177782],[-3,0,0.001757044],[-2.5,0,0.00254737],[-2,0,0.003589161],[-1.5,0,0.00491457],[-1,0,0.006539878],[-0.5,0,0.008457565],[0,0,0.0106295],[0.5,0,0.0129829],[1,0,0.01541069],[1.5,0,0.01777723],[2,0,0.01992956],[2.5,0,0.02171316],[3,0,0.02299005],[3.5,0,0.02365638],[4,0,0.02365638],[4.5,0,0.02299005],[5,0,0.02171316],[5.5,0,0.01992956],[6,0,0.01777723],[6.5,0,0.01541069],[7,0,0.0129829],[7.5,0,0.0106295],[8,0,0.008457565],[8.5,0,0.006539878],[9,0,0.00491457],[9.5,0,0.003589161],[10,0,0.00254737],[10.5,0,0.001757044],[11,0,0.001177782],[11.5,0,0.000767253],[12,0,0.0004857401],[12.5,0,0.0002988553],[13,0,0.0001786938],[13.5,0,0.0001038364],[14,0,5.86383e-005],[14.5,0,3.218138e-005],[15,0,1.716404e-005],[-5,0.5,0.000289402],[-4.5,0.5,0.0004737471],[-4,0.5,0.0007536737],[-3.5,0.5,0.00116523],[-3,0.5,0.00175078],[-2.5,0.5,0.002556484],[-2,0.5,0.003627823],[-1.5,0.5,0.005003118],[-1,0.5,0.006705435],[-0.5,0.5,0.008733831],[0,0.5,0.01105539],[0.5,0.5,0.01359988],[1,0.5,0.01625877],[1.5,0.5,0.01889],[2,0.5,0.02132886],[2.5,0.5,0.02340427],[3,0.5,0.02495824],[3.5,0.5,0.02586572],[4,0.5,0.02605113],[4.5,0.5,0.02549884],[5,0.5,0.02425524],[5.5,0.5,0.02242242],[6,0.5,0.02014424],[6.5,0.5,0.01758778],[7,0.5,0.01492322],[7.5,0.5,0.01230568],[8,0.5,0.009861445],[8.5,0.5,0.007680101],[9,0.5,0.005812794],[9.5,0.5,0.004275574],[10,0.5,0.003056297],[10.5,0.5,0.002123188],[11,0.5,0.001433417],[11.5,0.5,0.0009404778],[12,0.5,0.0005996751],[12.5,0.5,0.0003715995],[13,0.5,0.0002237823],[13.5,0.5,0.0001309688],[14,0.5,7.449066e-005],[14.5,0.5,4.11744e-005],[15,0.5,2.211792e-005],[-5,1,0.0002627986],[-4.5,1,0.0004332817],[-4,1,0.0006942392],[-3.5,1,0.001081035],[-3,1,0.001635918],[-2.5,1,0.002405887],[-2,1,0.003438589],[-1.5,1,0.00477614],[-1,1,0.006447115],[-0.5,1,0.008457565],[0,1,0.01078244],[0.5,1,0.01335919],[1,1,0.0160855],[1.5,1,0.01882266],[2,1,0.02140518],[2.5,1,0.02365638],[3,1,0.02540793],[3.5,1,0.02652051],[4,1,0.0269021],[4.5,1,0.02652051],[5,1,0.02540793],[5.5,1,0.02365638],[6,1,0.02140518],[6.5,1,0.01882266],[7,1,0.0160855],[7.5,1,0.01335919],[8,1,0.01078244],[8.5,1,0.008457565],[9,1,0.006447115],[9.5,1,0.00477614],[10,1,0.003438589],[10.5,1,0.002405887],[11,1,0.001635918],[11.5,1,0.001081035],[12,1,0.0006942392],[12.5,1,0.0004332817],[13,1,0.0002627986],[13.5,1,0.0001549057],[14,1,8.873671e-005],[14.5,1,4.940045e-005],[15,1,2.672699e-005],[-5,1.5,0.0002237823],[-4.5,1.5,0.0003715995],[-4,1.5,0.0005996751],[-3.5,1.5,0.0009404778],[-3,1.5,0.001433417],[-2.5,1.5,0.002123188],[-2,1.5,0.003056297],[-1.5,1.5,0.004275574],[-1,1.5,0.005812794],[-0.5,1.5,0.007680101],[0,1.5,0.009861445],[0.5,1.5,0.01230568],[1,1.5,0.01492322],[1.5,1.5,0.01758778],[2,1.5,0.02014424],[2.5,1.5,0.02242242],[3,1.5,0.02425524],[3.5,1.5,0.02549884],[4,1.5,0.02605113],[4.5,1.5,0.02586572],[5,1.5,0.02495824],[5.5,1.5,0.02340427],[6,1.5,0.02132886],[6.5,1.5,0.01889],[7,1.5,0.01625877],[7.5,1.5,0.01359988],[8,1.5,0.01105539],[8.5,1.5,0.008733831],[9,1.5,0.006705435],[9.5,1.5,0.005003118],[10,1.5,0.003627823],[10.5,1.5,0.002556484],[11,1.5,0.00175078],[11.5,1.5,0.00116523],[12,1.5,0.0007536737],[12.5,1.5,0.0004737471],[13,1.5,0.000289402],[13.5,1.5,0.0001718098],[14,1.5,9.912563e-005],[14.5,1.5,5.557962e-005],[15,1.5,3.028564e-005],[-5,2,0.0001786938],[-4.5,2,0.0002988553],[-4,2,0.0004857401],[-3.5,2,0.000767253],[-3,2,0.001177782],[-2.5,2,0.001757044],[-2,2,0.00254737],[-1.5,2,0.003589161],[-1,2,0.00491457],[-0.5,2,0.006539878],[0,2,0.008457565],[0.5,2,0.0106295],[1,2,0.0129829],[1.5,2,0.01541069],[2,2,0.01777723],[2.5,2,0.01992956],[3,2,0.02171316],[3.5,2,0.02299005],[4,2,0.02365638],[4.5,2,0.02365638],[5,2,0.02299005],[5.5,2,0.02171316],[6,2,0.01992956],[6.5,2,0.01777723],[7,2,0.01541069],[7.5,2,0.0129829],[8,2,0.0106295],[8.5,2,0.008457565],[9,2,0.006539878],[9.5,2,0.00491457],[10,2,0.003589161],[10.5,2,0.00254737],[11,2,0.001757044],[11.5,2,0.001177782],[12,2,0.000767253],[12.5,2,0.0004857401],[13,2,0.0002988553],[13.5,2,0.0001786938],[14,2,0.0001038364],[14.5,2,5.86383e-005],[15,2,3.218138e-005],[-5,2.5,0.0001338056],[-4.5,2.5,0.0002253865],[-4,2.5,0.0003689546],[-3.5,2.5,0.0005869616],[-3,2.5,0.000907482],[-2.5,2.5,0.001363509],[-2,2.5,0.001990991],[-1.5,2.5,0.002825351],[-1,2.5,0.003896432],[-0.5,2.5,0.005222199],[0,2.5,0.006801915],[0.5,2.5,0.00860995],[1,2.5,0.0105916],[1.5,2.5,0.01266235],[2,2.5,0.01471155],[2.5,2.5,0.01661094],[3,2.5,0.01822726],[3.5,2.5,0.0194375],[4,2.5,0.02014424],[4.5,2.5,0.02028864],[5,2.5,0.01985851],[5.5,2.5,0.01889],[6,2.5,0.0174626],[6.5,2.5,0.01568835],[7,2.5,0.01369737],[7.5,2.5,0.01162222],[8,2.5,0.009583677],[8.5,2.5,0.007680101],[9,2.5,0.005981269],[9.5,2.5,0.004527009],[10,2.5,0.003329821],[10.5,2.5,0.002380247],[11,2.5,0.00165354],[11.5,2.5,0.001116347],[12,2.5,0.0007324449],[12.5,2.5,0.0004670275],[13,2.5,0.000289402],[13.5,2.5,0.0001742818],[14,2.5,0.0001019986],[14.5,2.5,5.801339e-005],[15,2.5,3.206665e-005],[-5,3,9.395506e-005],[-4.5,3,0.0001593954],[-4,3,0.0002627986],[-3.5,3,0.0004210774],[-3,3,0.0006556806],[-2.5,3,0.0009922344],[-2,3,0.001459244],[-1.5,3,0.00208561],[-1,3,0.002896876],[-0.5,3,0.003910373],[0,3,0.005129773],[0.5,3,0.006539878],[1,3,0.008102755],[1.5,3,0.009756351],[2,3,0.01141652],[2.5,3,0.0129829],[3,3,0.01434832],[3.5,3,0.01541069],[4,3,0.0160855],[4.5,3,0.01631695],[5,3,0.0160855],[5.5,3,0.01541069],[6,3,0.01434832],[6.5,3,0.0129829],[7,3,0.01141652],[7.5,3,0.009756351],[8,3,0.008102755],[8.5,3,0.006539878],[9,3,0.005129773],[9.5,3,0.003910373],[10,3,0.002896876],[10.5,3,0.00208561],[11,3,0.001459244],[11.5,3,0.0009922344],[12,3,0.0006556806],[12.5,3,0.0004210774],[13,3,0.0002627986],[13.5,3,0.0001593954],[14,3,9.395506e-005],[14.5,3,5.382154e-005],[15,3,2.996288e-005],[-5,3.5,6.18653e-005],[-4.5,3.5,0.0001057073],[-4,3.5,0.0001755312],[-3.5,3.5,0.0002832665],[-3,3.5,0.0004442502],[-2.5,3.5,0.0006770984],[-2,3.5,0.001002923],[-1.5,3.5,0.001443693],[-1,3.5,0.002019638],[-0.5,3.5,0.002745769],[0,3.5,0.003627823],[0.5,3.5,0.004658217],[1,3.5,0.005812794],[1.5,3.5,0.00704923],[2,3.5,0.008307878],[2.5,3.5,0.009515466],[3,3.5,0.0105916],[3.5,3.5,0.01145737],[4,3.5,0.0120448],[4.5,3.5,0.01230568],[5,3.5,0.0122181],[5.5,3.5,0.01178944],[6,3.5,0.01105539],[6.5,3.5,0.01007504],[7,3.5,0.008923004],[7.5,3.5,0.007680101],[8,3.5,0.006424131],[8.5,3.5,0.005222199],[9,3.5,0.00412557],[9.5,3.5,0.003167423],[10,3.5,0.002363306],[10.5,3.5,0.001713662],[11,3.5,0.001207597],[11.5,3.5,0.0008270098],[12,3.5,0.0005504157],[12.5,3.5,0.0003560102],[13,3.5,0.0002237823],[13.5,3.5,0.0001367038],[14,3.5,8.11572e-005],[14.5,3.5,4.682363e-005],[15,3.5,2.625395e-005],[-5,4,3.819928e-005],[-4.5,4,6.573778e-005],[-4,4,0.0001099427],[-3.5,4,0.0001786938],[-3,4,0.0002822566],[-2.5,4,0.0004332817],[-2,4,0.0006463803],[-1.5,4,0.0009371249],[-1,4,0.001320379],[-0.5,4,0.001807969],[0,4,0.002405887],[0.5,4,0.003111364],[1,4,0.003910373],[1.5,4,0.00477614],[2,4,0.005669275],[2.5,4,0.006539878],[3,4,0.007331676],[3.5,4,0.007987824],[4,4,0.008457565],[4.5,4,0.008702695],[5,4,0.008702695],[5.5,4,0.008457565],[6,4,0.007987824],[6.5,4,0.007331676],[7,4,0.006539878],[7.5,4,0.005669275],[8,4,0.00477614],[8.5,4,0.003910373],[9,4,0.003111364],[9.5,4,0.002405887],[10,4,0.001807969],[10.5,4,0.001320379],[11,4,0.0009371249],[11.5,4,0.0006463803],[12,4,0.0004332817],[12.5,4,0.0002822566],[13,4,0.0001786938],[13.5,4,0.0001099427],[14,4,6.573778e-005],[14.5,4,3.819928e-005],[15,4,2.157182e-005],[-5,4.5,2.211792e-005],[-4.5,4.5,3.833595e-005],[-4,4.5,6.457431e-005],[-3.5,4.5,0.0001057073],[-3,4.5,0.0001681673],[-2.5,4.5,0.000259998],[-2,4.5,0.0003906518],[-1.5,4.5,0.0005704286],[-1,4.5,0.0008094767],[-0.5,4.5,0.001116347],[0,4.5,0.001496185],[0.5,4.5,0.001948781],[1,4.5,0.002466792],[1.5,4.5,0.003034544],[2,4.5,0.003627823],[2.5,4.5,0.004214929],[3,4.5,0.004759113],[3.5,4.5,0.005222199],[4,4.5,0.005568937],[4.5,4.5,0.005771422],[5,4.5,0.005812794],[5.5,4.5,0.005689559],[6,4.5,0.005412076],[6.5,4.5,0.005003118],[7,4.5,0.004494788],[7.5,4.5,0.003924364],[8,4.5,0.003329821],[8.5,4.5,0.002745769],[9,4.5,0.002200386],[9.5,4.5,0.001713662],[10,4.5,0.00129701],[10.5,4.5,0.0009540096],[11,4.5,0.0006819521],[11.5,4.5,0.0004737471],[12,4.5,0.0003198386],[12.5,4.5,0.000209849],[13,4.5,0.0001338056],[13.5,4.5,8.291505e-005],[14,4.5,4.993258e-005],[14.5,4.5,2.92231e-005],[15,4.5,1.662111e-005],[-5,5,1.200921e-005],[-4.5,5,2.096421e-005],[-4,5,3.556592e-005],[-3.5,5,5.86383e-005],[-3,5,9.395506e-005],[-2.5,5,0.0001463021],[-2,5,0.0002213974],[-1.5,5,0.0003256013],[-1,5,0.0004653625],[-0.5,5,0.0006463803],[0,5,0.0008725221],[0.5,5,0.001144607],[1,5,0.001459244],[1.5,5,0.001807969],[2,5,0.002176936],[2.5,5,0.00254737],[3,5,0.002896876],[3.5,5,0.003201543],[4,5,0.003438589],[4.5,5,0.003589161],[5,5,0.003640803],[5.5,5,0.003589161],[6,5,0.003438589],[6.5,5,0.003201543],[7,5,0.002896876],[7.5,5,0.00254737],[8,5,0.002176936],[8.5,5,0.001807969],[9,5,0.001459244],[9.5,5,0.001144607],[10,5,0.0008725221],[10.5,5,0.0006463803],[11,5,0.0004653625],[11.5,5,0.0003256013],[12,5,0.0002213974],[12.5,5,0.0001463021],[13,5,9.395506e-005],[13.5,5,5.86383e-005],[14,5,3.556592e-005],[14.5,5,2.096421e-005],[15,5,1.200921e-005],[-5,5.5,6.114566e-006],[-4.5,5.5,1.075058e-005],[-4,5.5,1.836917e-005],[-3.5,5.5,3.050274e-005],[-3,5.5,4.922433e-005],[-2.5,5.5,7.719912e-005],[-2,5.5,0.0001176621],[-1.5,5.5,0.0001742818],[-1,5.5,0.0002508762],[-0.5,5.5,0.0003509605],[0,5.5,0.0004771432],[0.5,5.5,0.0006304211],[1,5.5,0.0008094767],[1.5,5.5,0.001010112],[2,5.5,0.001224973],[2.5,5.5,0.001443693],[3,5.5,0.00165354],[3.5,5.5,0.001840544],[4,5.5,0.001990991],[4.5,5.5,0.002093072],[5,5.5,0.002138407],[5.5,5.5,0.002123188],[6,5.5,0.002048697],[6.5,5.5,0.001921139],[7,5.5,0.00175078],[7.5,5.5,0.001550586],[8,5.5,0.001334602],[8.5,5.5,0.001116347],[9,5.5,0.000907482],[9.5,5.5,0.0007169166],[10,5.5,0.0005504157],[10.5,5.5,0.000410681],[11,5.5,0.0002977898],[11.5,5.5,0.000209849],[12,5.5,0.0001437128],[12.5,5.5,9.56479e-005],[13,5.5,6.18653e-005],[13.5,5.5,3.888753e-005],[14,5.5,2.375556e-005],[14.5,5.5,1.410301e-005],[15,5.5,8.136727e-006],[-5,6,2.919429e-006],[-4.5,6,5.16971e-006],[-4,6,8.896641e-006],[-3.5,6,1.487913e-005],[-3,6,2.418358e-005],[-2.5,6,3.819928e-005],[-2,6,5.86383e-005],[-1.5,6,8.747806e-005],[-1,6,0.0001268261],[-0.5,6,0.0001786938],[0,6,0.000244682],[0.5,6,0.0003256013],[1,6,0.0004210774],[1.5,6,0.0005292114],[2,6,0.0006463803],[2.5,6,0.000767253],[3,6,0.0008850762],[3.5,6,0.0009922344],[4,6,0.001081035],[4.5,6,0.001144607],[5,6,0.001177782],[5.5,6,0.001177782],[6,6,0.001144607],[6.5,6,0.001081035],[7,6,0.0009922344],[7.5,6,0.0008850762],[8,6,0.000767253],[8.5,6,0.0006463803],[9,6,0.0005292114],[9.5,6,0.0004210774],[10,6,0.0003256013],[10.5,6,0.000244682],[11,6,0.0001786938],[11.5,6,0.0001268261],[12,6,8.747806e-005],[12.5,6,5.86383e-005],[13,6,3.819928e-005],[13.5,6,2.418358e-005],[14,6,1.487913e-005],[14.5,6,8.896641e-006],[15,6,5.16971e-006],[-5,6.5,1.307107e-006],[-4.5,6.5,2.331211e-006],[-4,6.5,4.040579e-006],[-3.5,6.5,6.806082e-006],[-3,6.5,1.114146e-005],[-2.5,6.5,1.772471e-005],[-2,6.5,2.740358e-005],[-1.5,6.5,4.11744e-005],[-1,6.5,6.012273e-005],[-0.5,6.5,8.531822e-005],[0,6.5,0.0001176621],[0.5,6.5,0.0001576967],[1,6.5,0.0002054],[1.5,6.5,0.000259998],[2,6.5,0.0003198386],[2.5,6.5,0.0003823697],[3,6.5,0.0004442502],[3.5,6.5,0.0005016068],[4,6.5,0.0005504157],[4.5,6.5,0.0005869616],[5,6.5,0.0006083034],[5.5,6.5,0.000612664],[6,6.5,0.0005996751],[6.5,6.5,0.0005704286],[7,6.5,0.0005273247],[7.5,6.5,0.0004737471],[8,6.5,0.0004136249],[8.5,6.5,0.0003509605],[9,6.5,0.000289402],[9.5,6.5,0.000231919],[10,6.5,0.0001806187],[10.5,6.5,0.0001367038],[11,6.5,0.0001005519],[11.5,6.5,7.187723e-005],[12,6.5,4.993258e-005],[12.5,6.5,3.371074e-005],[13,6.5,2.211792e-005],[13.5,6.5,1.410301e-005],[14,6.5,8.739182e-006],[14.5,6.5,5.262856e-006],[15,6.5,3.080092e-006],[-5,7,5.487894e-007],[-4.5,7,9.857758e-007],[-4,7,1.720847e-006],[-3.5,7,2.919429e-006],[-3,7,4.813325e-006],[-2.5,7,7.712301e-006],[-2,7,1.200921e-005],[-1.5,7,1.817341e-005],[-1,7,2.672699e-005],[-0.5,7,3.819928e-005],[0,7,5.305813e-005],[0.5,7,7.162098e-005],[1,7,9.395506e-005],[1.5,7,0.000119782],[2,7,0.0001484071],[2.5,7,0.0001786938],[3,7,0.0002091008],[3.5,7,0.00023779],[4,7,0.0002627986],[4.5,7,0.0002822566],[5,7,0.0002946163],[5.5,7,0.0002988553],[6,7,0.0002946163],[6.5,7,0.0002822566],[7,7,0.0002627986],[7.5,7,0.00023779],[8,7,0.0002091008],[8.5,7,0.0001786938],[9,7,0.0001484071],[9.5,7,0.000119782],[10,7,9.395506e-005],[10.5,7,7.162098e-005],[11,7,5.305813e-005],[11.5,7,3.819928e-005],[12,7,2.672699e-005],[12.5,7,1.817341e-005],[13,7,1.200921e-005],[13.5,7,7.712301e-006],[14,7,4.813325e-006],[14.5,7,2.919429e-006],[15,7,1.720847e-006],[-5,7.5,2.160634e-007],[-4.5,7.5,3.908911e-007],[-4,7.5,6.872614e-007],[-3.5,7.5,1.174302e-006],[-3,7.5,1.949975e-006],[-2.5,7.5,3.146806e-006],[-2,7.5,4.935174e-006],[-1.5,7.5,7.521884e-006],[-1,7.5,1.114146e-005],[-0.5,7.5,1.603798e-005],[0,7.5,2.243616e-005],[0.5,7.5,3.050274e-005],[1,7.5,4.030147e-005],[1.5,7.5,5.174812e-005],[2,7.5,6.457431e-005],[2.5,7.5,7.830988e-005],[3,7.5,9.229218e-005],[3.5,7.5,0.0001057073],[4,7.5,0.0001176621],[4.5,7.5,0.0001272798],[5,7.5,0.0001338056],[5.5,7.5,0.0001367038],[6,7.5,0.0001357308],[6.5,7.5,0.0001309688],[7,7.5,0.0001228143],[7.5,7.5,0.0001119236],[8,7.5,9.912563e-005],[8.5,7.5,8.531822e-005],[9,7.5,7.136565e-005],[9.5,7.5,5.801339e-005],[10,7.5,4.583094e-005],[10.5,7.5,3.518689e-005],[11,7.5,2.625395e-005],[11.5,7.5,1.903707e-005],[12,7.5,1.34152e-005],[12.5,7.5,9.18725e-006],[13,7.5,6.114566e-006],[13.5,7.5,3.954916e-006],[14,7.5,2.485997e-006],[14.5,7.5,1.518642e-006],[15,7.5,9.01575e-007],[-5,8,7.976966e-008],[-4.5,8,1.453498e-007],[-4,8,2.573847e-007],[-3.5,8,4.429376e-007],[-3,8,7.407882e-007],[-2.5,8,1.204029e-006],[-2,8,1.90183e-006],[-1.5,8,2.919429e-006],[-1,8,4.355276e-006],[-0.5,8,6.314298e-006],[0,8,8.896641e-006],[0.5,8,1.2182e-005],[1,8,1.621074e-005],[1.5,8,2.096421e-005],[2,8,2.634789e-005],[2.5,8,3.218138e-005],[3,8,3.819928e-005],[3.5,8,4.406535e-005],[4,8,4.940045e-005],[4.5,8,5.382154e-005],[5,8,5.698662e-005],[5.5,8,5.86383e-005],[6,8,5.86383e-005],[6.5,8,5.698662e-005],[7,8,5.382154e-005],[7.5,8,4.940045e-005],[8,8,4.406535e-005],[8.5,8,3.819928e-005],[9,8,3.218138e-005],[9.5,8,2.634789e-005],[10,8,2.096421e-005],[10.5,8,1.621074e-005],[11,8,1.2182e-005],[11.5,8,8.896641e-006],[12,8,6.314298e-006],[12.5,8,4.355276e-006],[13,8,2.919429e-006],[13.5,8,1.90183e-006],[14,8,1.204029e-006],[14.5,8,7.407882e-007],[15,8,4.429376e-007],[-5,8.5,2.761693e-008],[-4.5,8.5,5.068205e-008],[-4,8.5,9.039088e-008],[-3.5,8.5,1.566703e-007],[-3,8.5,2.639004e-007],[-2.5,8.5,4.320015e-007],[-2,8.5,6.872614e-007],[-1.5,8.5,1.062552e-006],[-1,8.5,1.596505e-006],[-0.5,8.5,2.331211e-006],[0,8.5,3.308146e-006],[0.5,8.5,4.562253e-006],[1,8.5,6.114566e-006],[1.5,8.5,7.964223e-006],[2,8.5,1.008121e-005],[2.5,8.5,1.240149e-005],[3,8.5,1.482608e-005],[3.5,8.5,1.722545e-005],[4,8.5,1.944941e-005],[4.5,8.5,2.134193e-005],[5,8.5,2.275897e-005],[5.5,8.5,2.358648e-005],[6,8.5,2.375556e-005],[6.5,8.5,2.325193e-005],[7,8.5,2.211792e-005],[7.5,8.5,2.04466e-005],[8,8.5,1.836917e-005],[8.5,8.5,1.603798e-005],[9,8.5,1.360822e-005],[9.5,8.5,1.122133e-005],[10,8.5,8.992474e-006],[10.5,8.5,7.003346e-006],[11,8.5,5.300582e-006],[11.5,8.5,3.898819e-006],[12,8.5,2.786983e-006],[12.5,8.5,1.936096e-006],[13,8.5,1.307107e-006],[13.5,8.5,8.576047e-007],[14,8.5,5.468329e-007],[14.5,8.5,3.388549e-007],[15,8.5,2.040631e-007],[-5,9,8.965904e-009],[-4.5,9,1.6572e-008],[-4,9,2.976785e-008],[-3.5,9,5.196508e-008],[-3,9,8.815911e-008],[-2.5,9,1.453498e-007],[-2,9,2.328913e-007],[-1.5,9,3.626467e-007],[-1,9,5.487894e-007],[-0.5,9,8.07085e-007],[0,9,1.153518e-006],[0.5,9,1.602217e-006],[1,9,2.162766e-006],[1.5,9,2.837197e-006],[2,9,3.617104e-006],[2.5,9,4.481507e-006],[3,9,5.396085e-006],[3.5,9,6.314298e-006],[4,9,7.180637e-006],[4.5,9,7.935831e-006],[5,9,8.523411e-006],[5.5,9,8.896641e-006],[6,9,9.024648e-006],[6.5,9,8.896641e-006],[7,9,8.523411e-006],[7.5,9,7.935831e-006],[8,9,7.180637e-006],[8.5,9,6.314298e-006],[9,9,5.396085e-006],[9.5,9,4.481507e-006],[10,9,3.617104e-006],[10.5,9,2.837197e-006],[11,9,2.162766e-006],[11.5,9,1.602217e-006],[12,9,1.153518e-006],[12.5,9,8.07085e-007],[13,9,5.487894e-007],[13.5,9,3.626467e-007],[14,9,2.328913e-007],[14.5,9,1.453498e-007],[15,9,8.815911e-008],[-5,9.5,2.729568e-009],[-4.5,9.5,5.081324e-009],[-4,9.5,9.192878e-009],[-3.5,9.5,1.616284e-008],[-3,9.5,2.761693e-008],[-2.5,9.5,4.585902e-008],[-2,9.5,7.400579e-008],[-1.5,9.5,1.160642e-007],[-1,9.5,1.768978e-007],[-0.5,9.5,2.620222e-007],[0,9.5,3.771771e-007],[0.5,9.5,5.276478e-007],[1,9.5,7.173558e-007],[1.5,9.5,9.477998e-007],[2,9.5,1.216999e-006],[2.5,9.5,1.518642e-006],[3,9.5,1.841672e-006],[3.5,9.5,2.170504e-006],[4,9.5,2.485997e-006],[4.5,9.5,2.767146e-006],[5,9.5,2.993335e-006],[5.5,9.5,3.146806e-006],[6,9.5,3.214965e-006],[6.5,9.5,3.192083e-006],[7,9.5,3.080092e-006],[7.5,9.5,2.888316e-006],[8,9.5,2.632191e-006],[8.5,9.5,2.331211e-006],[9,9.5,2.006492e-006],[9.5,9.5,1.678359e-006],[10,9.5,1.364344e-006],[10.5,9.5,1.07784e-006],[11,9.5,8.275165e-007],[11.5,9.5,6.174338e-007],[12,9.5,4.477089e-007],[12.5,9.5,3.154952e-007],[13,9.5,2.160634e-007],[13.5,9.5,1.438008e-007],[14,9.5,9.301072e-008],[14.5,9.5,5.846504e-008],[15,9.5,3.571504e-008],[-5,10,7.792463e-010],[-4.5,10,1.461032e-009],[-4,10,2.662175e-009],[-3.5,10,4.714166e-009],[-3,10,8.112687e-009],[-2.5,10,1.356801e-008],[-2,10,2.205257e-008],[-1.5,10,3.483323e-008],[-1,10,5.347121e-008],[-0.5,10,7.976966e-008],[0,10,1.156504e-007],[0.5,10,1.629477e-007],[1,10,2.231211e-007],[1.5,10,2.9691e-007],[2,10,3.839729e-007],[2.5,10,4.825783e-007],[3,10,5.894225e-007],[3.5,10,6.996441e-007],[4,10,8.07085e-007],[4.5,10,9.048007e-007],[5,10,9.857758e-007],[5.5,10,1.043746e-006],[6,10,1.073998e-006],[6.5,10,1.073998e-006],[7,10,1.043746e-006],[7.5,10,9.857758e-007],[8,10,9.048007e-007],[8.5,10,8.07085e-007],[9,10,6.996441e-007],[9.5,10,5.894225e-007],[10,10,4.825783e-007],[10.5,10,3.839729e-007],[11,10,2.9691e-007],[11.5,10,2.231211e-007],[12,10,1.629477e-007],[12.5,10,1.156504e-007],[13,10,7.976966e-008],[13.5,10,5.347121e-008],[14,10,3.483323e-008],[14.5,10,2.205257e-008],[15,10,1.356801e-008],[-5,10.5,2.086107e-010],[-4.5,10.5,3.939344e-010],[-4,10.5,7.229407e-010],[-3.5,10.5,1.289357e-009],[-3,10.5,2.234781e-009],[-2.5,10.5,3.764337e-009],[-2,10.5,6.16217e-009],[-1.5,10.5,9.803258e-009],[-1,10.5,1.515649e-008],[-0.5,10.5,2.277291e-008],[0,10.5,3.325295e-008],[0.5,10.5,4.718817e-008],[1,10.5,6.507705e-008],[1.5,10.5,8.72196e-008],[2,10.5,1.136035e-007],[2.5,10.5,1.438008e-007],[3,10.5,1.768978e-007],[3.5,10.5,2.114827e-007],[4,10.5,2.457079e-007],[4.5,10.5,2.774309e-007],[5,10.5,3.044263e-007],[5.5,10.5,3.246393e-007],[6,10.5,3.364431e-007],[6.5,10.5,3.388549e-007],[7,10.5,3.316709e-007],[7.5,10.5,3.154952e-007],[8,10.5,2.916551e-007],[8.5,10.5,2.620222e-007],[9,10.5,2.287695e-007],[9.5,10.5,1.941108e-007],[10,10.5,1.600637e-007],[10.5,10.5,1.282708e-007],[11,10.5,9.989736e-008],[11.5,10.5,7.560874e-008],[12,10.5,5.561367e-008],[12.5,10.5,3.975417e-008],[13,10.5,2.761693e-008],[13.5,10.5,1.864489e-008],[14,10.5,1.223307e-008],[14.5,10.5,7.800153e-009],[15,10.5,4.833505e-009],[-5,11,5.236965e-011],[-4.5,11,9.960221e-011],[-4,11,1.840983e-010],[-3.5,11,3.30691e-010],[-3,11,5.772798e-010],[-2.5,11,9.793593e-010],[-2,11,1.61469e-009],[-1.5,11,2.587189e-009],[-1,11,4.028641e-009],[-0.5,11,6.096499e-009],[0,11,8.965904e-009],[0.5,11,1.281443e-008],[1,11,1.779902e-008],[1.5,11,2.402616e-008],[2,11,3.151841e-008],[2.5,11,4.01824e-008],[3,11,4.978505e-008],[3.5,11,5.994509e-008],[4,11,7.014552e-008],[4.5,11,7.976966e-008],[5,11,8.815911e-008],[5.5,11,9.468654e-008],[6,11,9.883275e-008],[6.5,11,1.002548e-007],[7,11,9.883275e-008],[7.5,11,9.468654e-008],[8,11,8.815911e-008],[8.5,11,7.976966e-008],[9,11,7.014552e-008],[9.5,11,5.994509e-008],[10,11,4.978505e-008],[10.5,11,4.01824e-008],[11,11,3.151841e-008],[11.5,11,2.402616e-008],[12,11,1.779902e-008],[12.5,11,1.281443e-008],[13,11,8.965904e-009],[13.5,11,6.096499e-009],[14,11,4.028641e-009],[14.5,11,2.587189e-009],[15,11,1.61469e-009],[-5,11.5,1.232831e-011],[-4.5,11.5,2.361539e-011],[-4,11.5,4.396207e-011],[-3.5,11.5,7.953398e-011],[-3,11.5,1.39836e-010],[-2.5,11.5,2.389333e-010],[-2,11.5,3.967583e-010],[-1.5,11.5,6.402755e-010],[-1,11.5,1.004152e-009],[-0.5,11.5,1.530465e-009],[0,11.5,2.266936e-009],[0.5,11.5,3.263221e-009],[1,11.5,4.565048e-009],[1.5,11.5,6.206343e-009],[2,11.5,8.200075e-009],[2.5,11.5,1.05291e-008],[3,11.5,1.313883e-008],[3.5,11.5,1.593358e-008],[4,11.5,1.877854e-008],[4.5,11.5,2.150809e-008],[5,11.5,2.394051e-008],[5.5,11.5,2.589742e-008],[6,11.5,2.722521e-008],[6.5,11.5,2.78149e-008],[7,11.5,2.761693e-008],[7.5,11.5,2.664802e-008],[8,11.5,2.498883e-008],[8.5,11.5,2.277291e-008],[9,11.5,2.016893e-008],[9.5,11.5,1.735956e-008],[10,11.5,1.452065e-008],[10.5,11.5,1.180389e-008],[11,11.5,9.325147e-009],[11.5,11.5,7.15942e-009],[12,11.5,5.341849e-009],[12.5,11.5,3.873441e-009],[13,11.5,2.729568e-009],[13.5,11.5,1.869315e-009],[14,11.5,1.244121e-009],[14.5,11.5,8.047004e-010],[15,11.5,5.058217e-010],[-5,12,2.721503e-012],[-4.5,12,5.25052e-012],[-4,12,9.844359e-012],[-3.5,12,1.793759e-011],[-3,12,3.17638e-011],[-2.5,12,5.466285e-011],[-2,12,9.142052e-011],[-1.5,12,1.48589e-010],[-1,12,2.347045e-010],[-0.5,12,3.602862e-010],[0,12,5.374838e-010],[0.5,12,7.792463e-010],[1,12,1.097933e-009],[1.5,12,1.503378e-009],[2,12,2.000564e-009],[2.5,12,2.587189e-009],[3,12,3.251587e-009],[3.5,12,3.971498e-009],[4,12,4.714166e-009],[4.5,12,5.438096e-009],[5,12,6.096499e-009],[5.5,12,6.642106e-009],[6,12,7.032708e-009],[6.5,12,7.236541e-009],[7,12,7.236541e-009],[7.5,12,7.032708e-009],[8,12,6.642106e-009],[8.5,12,6.096499e-009],[9,12,5.438096e-009],[9.5,12,4.714166e-009],[10,12,3.971498e-009],[10.5,12,3.251587e-009],[11,12,2.587189e-009],[11.5,12,2.000564e-009],[12,12,1.503378e-009],[12.5,12,1.097933e-009],[13,12,7.792463e-010],[13.5,12,5.374838e-010],[14,12,3.602862e-010],[14.5,12,2.347045e-010],[15,12,1.48589e-010],[-5,12.5,5.633717e-013],[-4.5,12.5,1.094689e-012],[-4,12.5,2.067178e-012],[-3.5,12.5,3.793645e-012],[-3,12.5,6.765922e-012],[-2.5,12.5,1.172706e-011],[-2,12.5,1.975343e-011],[-1.5,12.5,3.23361e-011],[-1,12.5,5.144277e-011],[-0.5,12.5,7.953398e-011],[0,12.5,1.195013e-010],[0.5,12.5,1.744954e-010],[1,12.5,2.476208e-010],[1.5,12.5,3.41493e-010],[2,12.5,4.576864e-010],[2.5,12.5,5.961367e-010],[3,12.5,7.545973e-010],[3.5,12.5,9.282742e-010],[4,12.5,1.10976e-009],[4.5,12.5,1.289357e-009],[5,12.5,1.455824e-009],[5.5,12.5,1.597483e-009],[6,12.5,1.703551e-009],[6.5,12.5,1.765491e-009],[7,12.5,1.778147e-009],[7.5,12.5,1.740449e-009],[8,12.5,1.655567e-009],[8.5,12.5,1.530465e-009],[9,12.5,1.374966e-009],[9.5,12.5,1.200472e-009],[10,12.5,1.0186e-009],[10.5,12.5,8.399372e-010],[11,12.5,6.731031e-010],[11.5,12.5,5.242133e-010],[12,12.5,3.967583e-010],[12.5,12.5,2.918338e-010],[13,12.5,2.086107e-010],[13.5,12.5,1.449204e-010],[14,12.5,9.783938e-011],[14.5,12.5,6.419328e-011],[15,12.5,4.093145e-011],[-5,13,1.093609e-013],[-4.5,13,2.140228e-013],[-4,13,4.070515e-013],[-3.5,13,7.523678e-013],[-3,13,1.351458e-012],[-2.5,13,2.359211e-012],[-2,13,4.002418e-012],[-1.5,13,6.598871e-012],[-1,13,1.057325e-011],[-0.5,13,1.646413e-011],[0,13,2.4915e-011],[0.5,13,3.66416e-011],[1,13,5.236965e-011],[1.5,13,7.274052e-011],[2,13,9.818943e-011],[2.5,13,1.288085e-010],[3,13,1.642163e-010],[3.5,13,2.034601e-010],[4,13,2.449819e-010],[4.5,13,2.866687e-010],[5,13,3.260004e-010],[5.5,13,3.602862e-010],[6,13,3.869623e-010],[6.5,13,4.039069e-010],[7,13,4.097184e-010],[7.5,13,4.039069e-010],[8,13,3.869623e-010],[8.5,13,3.602862e-010],[9,13,3.260004e-010],[9.5,13,2.866687e-010],[10,13,2.449819e-010],[10.5,13,2.034601e-010],[11,13,1.642163e-010],[11.5,13,1.288085e-010],[12,13,9.818943e-011],[12.5,13,7.274052e-011],[13,13,5.236965e-011],[13.5,13,3.66416e-011],[14,13,2.4915e-011],[14.5,13,1.646413e-011],[15,13,1.057325e-011],[-5,13.5,1.990722e-014],[-4.5,13.5,3.923832e-014],[-4,13.5,7.516261e-014],[-3.5,13.5,1.399216e-013],[-3,13.5,2.531392e-013],[-2.5,13.5,4.450672e-013],[-2,13.5,7.604722e-013],[-1.5,13.5,1.262795e-012],[-1,13.5,2.037857e-012],[-0.5,13.5,3.195995e-012],[0,13.5,4.871136e-012],[0.5,13.5,7.215159e-012],[1,13.5,1.038612e-011],[1.5,13.5,1.452955e-011],[2,13.5,1.975343e-011],[2.5,13.5,2.609904e-011],[3,13.5,3.351183e-011],[3.5,13.5,4.181801e-011],[4,13.5,5.07131e-011],[4.5,13.5,5.976798e-011],[5,13.5,6.845553e-011],[5.5,13.5,7.619739e-011],[6,13.5,8.242581e-011],[6.5,13.5,8.665187e-011],[7,13.5,8.852873e-011],[7.5,13.5,8.789864e-011],[8,13.5,8.481479e-011],[8.5,13.5,7.953398e-011],[9,13.5,7.248119e-011],[9.5,13.5,6.419328e-011],[10,13.5,5.525167e-011],[10.5,13.5,4.621605e-011],[11,13.5,3.756919e-011],[11.5,13.5,2.96799e-011],[12,13.5,2.278687e-011],[12.5,13.5,1.700193e-011],[13,13.5,1.232831e-011],[13.5,13.5,8.687617e-012],[14,13.5,5.949619e-012],[14.5,13.5,3.959764e-012],[15,13.5,2.561185e-012],[-5,14,3.39813e-015],[-4.5,14,6.745931e-015],[-4,14,1.301474e-014],[-3.5,14,2.440173e-014],[-3,14,4.446284e-014],[-2.5,14,7.873458e-014],[-2,14,1.354957e-013],[-1.5,14,2.266088e-013],[-1,14,3.683154e-013],[-0.5,14,5.817742e-013],[0,14,8.930601e-013],[0.5,14,1.332289e-012],[1,14,1.931558e-012],[1.5,14,2.721503e-012],[2,14,3.726502e-012],[2.5,14,4.958902e-012],[3,14,6.413e-012],[3.5,14,8.059879e-012],[4,14,9.844359e-012],[4.5,14,1.168525e-011],[5,14,1.347969e-011],[5.5,14,1.511171e-011],[6,14,1.646413e-011],[6.5,14,1.743234e-011],[7,14,1.793759e-011],[7.5,14,1.793759e-011],[8,14,1.743234e-011],[8.5,14,1.646413e-011],[9,14,1.511171e-011],[9.5,14,1.347969e-011],[10,14,1.168525e-011],[10.5,14,9.844359e-012],[11,14,8.059879e-012],[11.5,14,6.413e-012],[12,14,4.958902e-012],[12.5,14,3.726502e-012],[13,14,2.721503e-012],[13.5,14,1.931558e-012],[14,14,1.332289e-012],[14.5,14,8.930601e-013],[15,14,5.817742e-013],[-5,14.5,5.439393e-016],[-4.5,14.5,1.087563e-015],[-4,14.5,2.113246e-015],[-3.5,14.5,3.990592e-015],[-3,14.5,7.323457e-015],[-2.5,14.5,1.30613e-014],[-2,14.5,2.263854e-014],[-1.5,14.5,3.813309e-014],[-1,14.5,6.242336e-014],[-0.5,14.5,9.930792e-014],[0,14.5,1.535367e-013],[0.5,14.5,2.306918e-013],[1,14.5,3.368554e-013],[1.5,14.5,4.780206e-013],[2,14.5,6.592366e-013],[2.5,14.5,8.835426e-013],[3,14.5,1.150815e-012],[3.5,14.5,1.456716e-012],[4,14.5,1.791991e-012],[4.5,14.5,2.14234e-012],[5,14.5,2.489044e-012],[5.5,14.5,2.810401e-012],[6,14.5,3.083867e-012],[6.5,14.5,3.288627e-012],[7,14.5,3.4082e-012],[7.5,14.5,3.432632e-012],[8,14.5,3.359858e-012],[8.5,14.5,3.195995e-012],[9,14.5,2.954493e-012],[9.5,14.5,2.654309e-012],[10,14.5,2.317456e-012],[10.5,14.5,1.96636e-012],[11,14.5,1.62146e-012],[11.5,14.5,1.299395e-012],[12,14.5,1.01197e-012],[12.5,14.5,7.659236e-013],[13,14.5,5.633717e-013],[13.5,14.5,4.027135e-013],[14,14.5,2.797621e-013],[14.5,14.5,1.888744e-013],[15,14.5,1.239222e-013],[-5,15,8.164733e-017],[-4.5,15,1.644175e-016],[-4,15,3.217702e-016],[-3.5,15,6.119771e-016],[-3,15,1.131139e-015],[-2.5,15,2.031836e-015],[-2,15,3.54693e-015],[-1.5,15,6.017392e-015],[-1,15,9.921001e-015],[-0.5,15,1.589624e-014],[0,15,2.475282e-014],[0.5,15,3.745819e-014],[1,15,5.508842e-014],[1.5,15,7.873458e-014],[2,15,1.093609e-013],[2.5,15,1.476218e-013],[3,15,1.936558e-013],[3.5,15,2.468892e-013],[4,15,3.058899e-013],[4.5,15,3.683154e-013],[5,15,4.30989e-013],[5.5,15,4.901218e-013],[6,15,5.416683e-013],[6.5,15,5.817742e-013],[7,15,6.072494e-013],[7.5,15,6.159866e-013],[8,15,6.072494e-013],[8.5,15,5.817742e-013],[9,15,5.416683e-013],[9.5,15,4.901218e-013],[10,15,4.30989e-013],[10.5,15,3.683154e-013],[11,15,3.058899e-013],[11.5,15,2.468892e-013],[12,15,1.936558e-013],[12.5,15,1.476218e-013],[13,15,1.093609e-013],[13.5,15,7.873458e-014],[14,15,5.508842e-014],[14.5,15,3.745819e-014],[15,15,2.475282e-014]],"ignoreExtent":false,"flags":128},"171":{"id":171,"type":"text","material":{"lit":false},"vertices":[[5,-8.39,-0.004559905]],"colors":[[0,0,0,1]],"texts":[["Y1"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[5,-8.39,-0.004559905]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":4136},"172":{"id":172,"type":"text","material":{"lit":false},"vertices":[[-8.39,5,-0.004559905]],"colors":[[0,0,0,1]],"texts":[["Y2"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-8.39,5,-0.004559905]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":4136},"173":{"id":173,"type":"text","material":{"lit":false},"vertices":[[-8.39,-8.39,0.01345105]],"colors":[[0,0,0,1]],"texts":[["Density"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-8.39,-8.39,0.01345105]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":4136},"174":{"id":174,"type":"linestrip","material":{"lit":false},"vertices":[[-5,-5,1.146066e-005],[-4.5,-5,1.79456e-005],[-4,-5,2.826249e-005],[-3.5,-5,4.302606e-005],[-3,-5,6.224708e-005],[-2.5,-5,8.488809e-005],[-2,-5,0.0001095062],[-1.5,-5,0.0001354254],[-1,-5,0.0001627141],[-0.5,-5,0.0001904372],[0,-5,0.0002158209],[0.5,-5,0.0002363462],[1,-5,0.0002515781],[1.5,-5,0.0002620281],[2,-5,0.0002678616],[2.5,-5,0.0002711958],[3,-5,0.0002817815],[3.5,-5,0.0002963529],[4,-5,0.0002775255],[4.5,-5,0.0002341229],[5,-5,0.0002044563],[5.5,-5,0.0001868718],[6,-5,0.0001718011],[6.5,-5,0.0001562788],[7,-5,0.000139487],[7.5,-5,0.0001214974],[8,-5,0.0001024068],[8.5,-5,8.197113e-005],[9,-5,6.098348e-005],[9.5,-5,4.201412e-005],[10,-5,2.748056e-005],[10.5,-5,1.774059e-005],[11,-5,1.156565e-005],[11.5,-5,7.565056e-006],[12,-5,4.852432e-006],[12.5,-5,2.98488e-006],[13,-5,1.735515e-006],[13.5,-5,9.481955e-007],[14,-5,4.875333e-007],[14.5,-5,2.376798e-007],[15,-5,1.111036e-007],[-5,-4.5,2.237209e-005],[-4.5,-4.5,3.180327e-005],[-4,-4.5,4.912843e-005],[-3.5,-4.5,7.398456e-005],[-3,-4.5,0.0001062965],[-2.5,-4.5,0.0001444892],[-2,-4.5,0.0001864878],[-1.5,-4.5,0.0002318203],[-1,-4.5,0.0002811916],[-0.5,-4.5,0.0003324109],[0,-4.5,0.0003786405],[0.5,-4.5,0.0004147704],[1,-4.5,0.0004414355],[1.5,-4.5,0.0004599946],[2,-4.5,0.0004701946],[2.5,-4.5,0.0004742752],[3,-4.5,0.000483124],[3.5,-4.5,0.0004937478],[4,-4.5,0.0004665218],[4.5,-4.5,0.0004111622],[5,-4.5,0.0003687843],[5.5,-4.5,0.0003377087],[6,-4.5,0.0003081178],[6.5,-4.5,0.0002773233],[7,-4.5,0.0002449829],[7.5,-4.5,0.0002116368],[8,-4.5,0.0001777556],[8.5,-4.5,0.0001425613],[9,-4.5,0.0001065677],[9.5,-4.5,7.380409e-005],[10,-4.5,4.854332e-005],[10.5,-4.5,3.157794e-005],[11,-4.5,2.08178e-005],[11.5,-4.5,1.380867e-005],[12,-4.5,8.987108e-006],[12.5,-4.5,5.601876e-006],[13,-4.5,3.293227e-006],[13.5,-4.5,1.814594e-006],[14,-4.5,9.383666e-007],[14.5,-4.5,4.586728e-007],[15,-4.5,2.142398e-007],[-5,-4,4.112635e-005],[-4.5,-4,5.3383e-005],[-4,-4,8.109483e-005],[-3.5,-4,0.0001204543],[-3,-4,0.0001714598],[-2.5,-4,0.0002326409],[-2,-4,0.0003012976],[-1.5,-4,0.000377329],[-1,-4,0.0004629512],[-0.5,-4,0.0005543241],[0,-4,0.0006346446],[0.5,-4,0.0006932709],[1,-4,0.000737445],[1.5,-4,0.0007713854],[2,-4,0.0007912981],[2.5,-4,0.0007962267],[3,-4,0.0007920349],[3.5,-4,0.0007798162],[4,-4,0.0007454995],[4.5,-4,0.0006956958],[5,-4,0.0006462956],[5.5,-4,0.0005930346],[6,-4,0.0005372703],[6.5,-4,0.0004808864],[7,-4,0.0004220904],[7.5,-4,0.0003589195],[8,-4,0.0002956315],[8.5,-4,0.0002345179],[9,-4,0.00017526],[9.5,-4,0.0001223395],[10,-4,8.160929e-005],[10.5,-4,5.40173e-005],[11,-4,3.621635e-005],[11.5,-4,2.437729e-005],[12,-4,1.606446e-005],[12.5,-4,1.012191e-005],[13,-4,6.007861e-006],[13.5,-4,3.339113e-006],[14,-4,1.739856e-006],[14.5,-4,8.55642e-007],[15,-4,4.012541e-007],[-5,-3.5,6.327545e-005],[-4.5,-3.5,8.672065e-005],[-4,-3.5,0.000132443],[-3.5,-3.5,0.0001909186],[-3,-3.5,0.0002642531],[-2.5,-3.5,0.0003558708],[-2,-3.5,0.0004623005],[-1.5,-3.5,0.0005829262],[-1,-3.5,0.000723041],[-0.5,-3.5,0.0008791749],[0,-3.5,0.001015789],[0.5,-3.5,0.001107298],[1,-3.5,0.001176776],[1.5,-3.5,0.001239609],[2,-3.5,0.001280798],[2.5,-3.5,0.001292596],[3,-3.5,0.00128148],[3.5,-3.5,0.001253596],[4,-3.5,0.001213399],[4.5,-3.5,0.001173225],[5,-3.5,0.001118858],[5.5,-3.5,0.001016316],[6,-3.5,0.000907936],[6.5,-3.5,0.0008119827],[7,-3.5,0.0007128508],[7.5,-3.5,0.0005953113],[8,-3.5,0.0004763827],[8.5,-3.5,0.0003724136],[9,-3.5,0.0002785494],[9.5,-3.5,0.0001954714],[10,-3.5,0.0001319688],[10.5,-3.5,8.922184e-005],[11,-3.5,6.112921e-005],[11.5,-3.5,4.182365e-005],[12,-3.5,2.788037e-005],[12.5,-3.5,1.772779e-005],[13,-3.5,1.061292e-005],[13.5,-3.5,5.950761e-006],[14,-3.5,3.129319e-006],[14.5,-3.5,1.553438e-006],[15,-3.5,7.349277e-007],[-5,-3,9.086989e-005],[-4.5,-3,0.0001408084],[-4,-3,0.0002162213],[-3.5,-3,0.0002972678],[-3,-3,0.0003921731],[-2.5,-3,0.0005203511],[-2,-3,0.0006772354],[-1.5,-3,0.0008578757],[-1,-3,0.001066805],[-0.5,-3,0.001303766],[0,-3,0.001525707],[0.5,-3,0.001685078],[1,-3,0.001809782],[1.5,-3,0.001927267],[2,-3,0.002012582],[2.5,-3,0.002050395],[3,-3,0.002030791],[3.5,-3,0.001995545],[4,-3,0.001950376],[4.5,-3,0.00191211],[5,-3,0.001838514],[5.5,-3,0.001657118],[6,-3,0.001466172],[6.5,-3,0.001304361],[7,-3,0.001145016],[7.5,-3,0.0009816453],[8,-3,0.0007547182],[8.5,-3,0.0005843053],[9,-3,0.0004421417],[9.5,-3,0.0003100148],[10,-3,0.0002091989],[10.5,-3,0.0001438895],[11,-3,0.0001008534],[11.5,-3,7.020459e-005],[12,-3,4.734426e-005],[12.5,-3,3.037661e-005],[13,-3,1.834266e-005],[13.5,-3,1.037853e-005],[14,-3,5.51339e-006],[14.5,-3,2.769239e-006],[15,-3,1.327097e-006],[-5,-2.5,0.0001442305],[-4.5,-2.5,0.0002205195],[-4,-2.5,0.0003183437],[-3.5,-2.5,0.0004252557],[-3,-2.5,0.0005536033],[-2.5,-2.5,0.0007333626],[-2,-2.5,0.0009581725],[-1.5,-2.5,0.001216761],[-1,-2.5,0.001502542],[-0.5,-2.5,0.001815174],[0,-2.5,0.002135989],[0.5,-2.5,0.002436238],[1,-2.5,0.002771712],[1.5,-2.5,0.002912784],[2,-2.5,0.003119577],[2.5,-2.5,0.003412267],[3,-2.5,0.003185678],[3.5,-2.5,0.003097689],[4,-2.5,0.003042762],[4.5,-2.5,0.002960254],[5,-2.5,0.002793584],[5.5,-2.5,0.002527546],[6,-2.5,0.002248884],[6.5,-2.5,0.001980725],[7,-2.5,0.001736836],[7.5,-2.5,0.001621885],[8,-2.5,0.001175007],[8.5,-2.5,0.0008900007],[9,-2.5,0.0006791472],[9.5,-2.5,0.0004778789],[10,-2.5,0.0003248813],[10.5,-2.5,0.0002277697],[11,-2.5,0.0001631131],[11.5,-2.5,0.0001154858],[12,-2.5,7.891313e-005],[12.5,-2.5,5.122688e-005],[13,-2.5,3.126673e-005],[13.5,-2.5,1.786575e-005],[14,-2.5,9.590415e-006],[14.5,-2.5,4.880421e-006],[15,-2.5,2.376079e-006],[-5,-2,0.0002165681],[-4.5,-2,0.0003067722],[-4,-2,0.0004096671],[-3.5,-2,0.0005531394],[-3,-2,0.0007481863],[-2.5,-2,0.001009775],[-2,-2,0.001321894],[-1.5,-2,0.001676139],[-1,-2,0.002057783],[-0.5,-2,0.002455546],[0,-2,0.002884016],[0.5,-2,0.003378617],[1,-2,0.00399484],[1.5,-2,0.004209686],[2,-2,0.004477599],[2.5,-2,0.004908534],[3,-2,0.004767854],[3.5,-2,0.004671034],[4,-2,0.004610468],[4.5,-2,0.004455058],[5,-2,0.00412845],[5.5,-2,0.003718831],[6,-2,0.003315366],[6.5,-2,0.002901602],[7,-2,0.002481993],[7.5,-2,0.002087073],[8,-2,0.001655077],[8.5,-2,0.001289009],[9,-2,0.0009762045],[9.5,-2,0.0006985588],[10,-2,0.0004903682],[10.5,-2,0.0003524965],[11,-2,0.000256847],[11.5,-2,0.0001845717],[12,-2,0.0001280628],[12.5,-2,8.446218e-005],[13,-2,5.228769e-005],[13.5,-2,3.0224e-005],[14,-2,1.641317e-005],[14.5,-2,8.478482e-006],[15,-2,4.205202e-006],[-5,-1.5,0.0002572473],[-4.5,-1.5,0.000368599],[-4,-1.5,0.0005108458],[-3.5,-1.5,0.0007113838],[-3,-1.5,0.00098646],[-2.5,-1.5,0.001347442],[-2,-1.5,0.001758092],[-1.5,-1.5,0.002231156],[-1,-1.5,0.002761184],[-0.5,-1.5,0.003279756],[0,-1.5,0.003830162],[0.5,-1.5,0.00450843],[1,-1.5,0.005235973],[1.5,-1.5,0.005755024],[2,-1.5,0.00614239],[2.5,-1.5,0.006588121],[3,-1.5,0.00670407],[3.5,-1.5,0.006766361],[4,-1.5,0.006806706],[4.5,-1.5,0.006441725],[5,-1.5,0.005880353],[5.5,-1.5,0.005301803],[6,-1.5,0.00472252],[6.5,-1.5,0.004121446],[7,-1.5,0.003510275],[7.5,-1.5,0.002946882],[8,-1.5,0.002472274],[8.5,-1.5,0.001842721],[9,-1.5,0.001380014],[9.5,-1.5,0.001006347],[10,-1.5,0.0007270694],[10.5,-1.5,0.0005303404],[11,-1.5,0.000389386],[11.5,-1.5,0.0002825409],[12,-1.5,0.0001988005],[12.5,-1.5,0.0001332903],[13,-1.5,8.380409e-005],[13.5,-1.5,4.912484e-005],[14,-1.5,2.7098e-005],[14.5,-1.5,1.428011e-005],[15,-1.5,7.2475e-006],[-5,-1,0.0002865055],[-4.5,-1,0.0004298274],[-4,-1,0.0006351962],[-3.5,-1,0.0008965567],[-3,-1,0.001237694],[-2.5,-1,0.00169112],[-2,-1,0.002216013],[-1.5,-1,0.002850228],[-1,-1,0.003595983],[-0.5,-1,0.004278233],[0,-1,0.005005321],[0.5,-1,0.005856466],[1,-1,0.006774438],[1.5,-1,0.00756753],[2,-1,0.008159697],[2.5,-1,0.008697816],[3,-1,0.008972279],[3.5,-1,0.009346751],[4,-1,0.009260431],[4.5,-1,0.008641054],[5,-1,0.008031354],[5.5,-1,0.007312033],[6,-1,0.006482272],[6.5,-1,0.005621665],[7,-1,0.004788855],[7.5,-1,0.004025815],[8,-1,0.003372102],[8.5,-1,0.002587227],[9,-1,0.001964825],[9.5,-1,0.001459167],[10,-1,0.001069749],[10.5,-1,0.0007746161],[11,-1,0.0005653928],[11.5,-1,0.0004110206],[12,-1,0.0002910183],[12.5,-1,0.0001969081],[13,-1,0.0001252403],[13.5,-1,7.457702e-005],[14,-1,4.205952e-005],[14.5,-1,2.278823e-005],[15,-1,1.19076e-005],[-5,-0.5,0.0003415506],[-4.5,-0.5,0.0005112388],[-4,-0.5,0.0007535119],[-3.5,-0.5,0.001067786],[-3,-0.5,0.001480715],[-2.5,-0.5,0.00202831],[-2,-0.5,0.002674827],[-1.5,-0.5,0.003443498],[-1,-0.5,0.004380933],[-0.5,-0.5,0.005345846],[0,-0.5,0.006522307],[0.5,-0.5,0.00752635],[1,-0.5,0.008530652],[1.5,-0.5,0.009576255],[2,-0.5,0.01040816],[2.5,-0.5,0.01100683],[3,-0.5,0.01156096],[3.5,-0.5,0.01200045],[4,-0.5,0.01195963],[4.5,-0.5,0.01119982],[5,-0.5,0.0105005],[5.5,-0.5,0.009608634],[6,-0.5,0.0084919],[6.5,-0.5,0.007322086],[7,-0.5,0.006249747],[7.5,-0.5,0.005257108],[8,-0.5,0.004327205],[8.5,-0.5,0.003470191],[9,-0.5,0.002690308],[9.5,-0.5,0.002049894],[10,-0.5,0.001525061],[10.5,-0.5,0.001089799],[11,-0.5,0.0007881323],[11.5,-0.5,0.000570958],[12,-0.5,0.0004018823],[12.5,-0.5,0.0002705456],[13,-0.5,0.0001724871],[13.5,-0.5,0.0001042258],[14,-0.5,6.040592e-005],[14.5,-0.5,3.387763e-005],[15,-0.5,1.833233e-005],[-5,0,0.000404885],[-4.5,0,0.0006008795],[-4,0,0.0008751656],[-3.5,0,0.001247549],[-3,0,0.001736728],[-2.5,0,0.002367277],[-2,0,0.003155175],[-1.5,0,0.004007886],[-1,0,0.005100699],[-0.5,0,0.006445236],[0,0,0.007772099],[0.5,0,0.009231257],[1,0,0.0103401],[1.5,0,0.01150674],[2,0,0.01267288],[2.5,0,0.01361682],[3,0,0.01442646],[3.5,0,0.01490689],[4,0,0.01484609],[4.5,0,0.0139666],[5,0,0.01310446],[5.5,0,0.01197915],[6,0,0.01063787],[6.5,0,0.009245371],[7,0,0.007866147],[7.5,0,0.006665982],[8,0,0.005502604],[8.5,0,0.00439601],[9,0,0.003433482],[9.5,0,0.002644341],[10,0,0.001983844],[10.5,0,0.001441925],[11,0,0.001056572],[11.5,0,0.0007655464],[12,0,0.0005293552],[12.5,0,0.0003487939],[13,0,0.0002205494],[13.5,0,0.0001347026],[14,0,8.026281e-005],[14.5,0,4.670894e-005],[15,0,2.629063e-005],[-5,0.5,0.0004554605],[-4.5,0.5,0.0006792492],[-4,0.5,0.000992506],[-3.5,0.5,0.001417752],[-3,0.5,0.001964358],[-2.5,0.5,0.002657247],[-2,0.5,0.00366049],[-1.5,0.5,0.004595141],[-1,0.5,0.00576894],[-0.5,0.5,0.007459366],[0,0.5,0.00887116],[0.5,0.5,0.01045944],[1,0.5,0.01193677],[1.5,0.5,0.01327685],[2,0.5,0.01474313],[2.5,0.5,0.01608566],[3,0.5,0.0168963],[3.5,0.5,0.01734479],[4,0.5,0.01754291],[4.5,0.5,0.01664309],[5,0.5,0.01560924],[5.5,0.5,0.01441156],[6,0.5,0.0129435],[6.5,0.5,0.0114458],[7,0.5,0.00960853],[7.5,0.5,0.008114891],[8,0.5,0.006699379],[8.5,0.5,0.005356373],[9,0.5,0.004211412],[9.5,0.5,0.003240301],[10,0.5,0.002437766],[10.5,0.5,0.001810906],[11,0.5,0.0013517],[11.5,0.5,0.0009828239],[12,0.5,0.0006644871],[12.5,0.5,0.0004252456],[13,0.5,0.0002654321],[13.5,0.5,0.0001632694],[14,0.5,9.966295e-005],[14.5,0.5,6.011502e-005],[15,0.5,3.534319e-005],[-5,1,0.0004830987],[-4.5,1,0.0007234502],[-4,1,0.001059132],[-3.5,1,0.001505477],[-3,1,0.002079298],[-2.5,1,0.002811247],[-2,1,0.003819107],[-1.5,1,0.00497342],[-1,1,0.006137304],[-0.5,1,0.007744041],[0,1,0.009407871],[0.5,1,0.01132513],[1,1,0.01336012],[1.5,1,0.01477714],[2,1,0.01638014],[2.5,1,0.01806297],[3,1,0.01928366],[3.5,1,0.01978466],[4,1,0.01963364],[4.5,1,0.01880674],[5,1,0.01787322],[5.5,1,0.01652871],[6,1,0.01492306],[6.5,1,0.0134116],[7,1,0.01134544],[7.5,1,0.00951534],[8,1,0.007900327],[8.5,1,0.00638144],[9,1,0.005065174],[9.5,1,0.003906359],[10,1,0.002961425],[10.5,1,0.002188874],[11,1,0.001609593],[11.5,1,0.001159115],[12,1,0.0007749737],[12.5,1,0.0004902227],[13,1,0.0003044707],[13.5,1,0.0001883835],[14,1,0.0001171817],[14.5,1,7.308921e-005],[15,1,4.508184e-005],[-5,1.5,0.0004881235],[-4.5,1.5,0.0007284761],[-4,1.5,0.001065357],[-3.5,1.5,0.001517409],[-3,1.5,0.002097141],[-2.5,1.5,0.002801118],[-2,1.5,0.003797894],[-1.5,1.5,0.004950379],[-1,1.5,0.006141036],[-0.5,1.5,0.007881681],[0,1.5,0.009610836],[0.5,1.5,0.01153261],[1,1.5,0.01378538],[1.5,1.5,0.01543789],[2,1.5,0.01713712],[2.5,1.5,0.01884047],[3,1.5,0.02055656],[3.5,1.5,0.02144197],[4,1.5,0.02136368],[4.5,1.5,0.02052283],[5,1.5,0.01955397],[5.5,1.5,0.01827986],[6,1.5,0.0165303],[6.5,1.5,0.0145778],[7,1.5,0.01255322],[7.5,1.5,0.01068448],[8,1.5,0.00909933],[8.5,1.5,0.007381013],[9,1.5,0.005875655],[9.5,1.5,0.004506598],[10,1.5,0.00346514],[10.5,1.5,0.00251358],[11,1.5,0.001780023],[11.5,1.5,0.001254741],[12,1.5,0.0008444163],[12.5,1.5,0.000541629],[13,1.5,0.0003381055],[13.5,1.5,0.0002097826],[14,1.5,0.0001319739],[14.5,1.5,8.473112e-005],[15,1.5,5.482933e-005],[-5,2,0.0004730257],[-4.5,2,0.000700683],[-4,2,0.001024601],[-3.5,2,0.001507917],[-3,2,0.002236372],[-2.5,2,0.002774193],[-2,2,0.003647106],[-1.5,2,0.004689286],[-1,2,0.00595497],[-0.5,2,0.007754375],[0,2,0.009409109],[0.5,2,0.01114514],[1,2,0.01317165],[1.5,2,0.01518674],[2,2,0.01712975],[2.5,2,0.01882221],[3,2,0.0208555],[3.5,2,0.02175393],[4,2,0.02216771],[4.5,2,0.021446],[5,2,0.02035583],[5.5,2,0.01902198],[6,2,0.01742218],[6.5,2,0.01507541],[7,2,0.01309105],[7.5,2,0.01135501],[8,2,0.009636084],[8.5,2,0.007805859],[9,2,0.006227697],[9.5,2,0.005007782],[10,2,0.003783584],[10.5,2,0.002677933],[11,2,0.001877535],[11.5,2,0.001313652],[12,2,0.0009009013],[12.5,2,0.0005906412],[13,2,0.0003692903],[13.5,2,0.0002274542],[14,2,0.0001430757],[14.5,2,9.355474e-005],[15,2,6.283284e-005],[-5,2.5,0.000436783],[-4.5,2.5,0.0006396043],[-4,2.5,0.0009291347],[-3.5,2.5,0.001357811],[-3,2.5,0.001960082],[-2.5,2.5,0.002600539],[-2,2.5,0.00344066],[-1.5,2.5,0.004441527],[-1,2.5,0.005637107],[-0.5,2.5,0.007220118],[0,2.5,0.008836407],[0.5,2.5,0.01041832],[1,2.5,0.01233387],[1.5,2.5,0.01431841],[2,2.5,0.01620566],[2.5,2.5,0.01812897],[3,2.5,0.01996089],[3.5,2.5,0.02129149],[4,2.5,0.02176238],[4.5,2.5,0.02133545],[5,2.5,0.02015966],[5.5,2.5,0.01897212],[6,2.5,0.0177619],[6.5,2.5,0.01513214],[7,2.5,0.01288958],[7.5,2.5,0.0112167],[8,2.5,0.009383049],[8.5,2.5,0.007482089],[9,2.5,0.006160055],[9.5,2.5,0.005133534],[10,2.5,0.003891281],[10.5,2.5,0.002751673],[11,2.5,0.001946777],[11.5,2.5,0.00136052],[12,2.5,0.0009550565],[12.5,2.5,0.0006391706],[13,2.5,0.0003964631],[13.5,2.5,0.0002396129],[14,2.5,0.0001488367],[14.5,2.5,9.741767e-005],[15,2.5,6.643661e-005],[-5,3,0.0004756195],[-4.5,3,0.0005826598],[-4,3,0.0007925265],[-3.5,3,0.001162311],[-3,3,0.001679428],[-2.5,3,0.002329326],[-2,3,0.003104219],[-1.5,3,0.004049969],[-1,3,0.00525702],[-0.5,3,0.006698369],[0,3,0.007916006],[0.5,3,0.009484243],[1,3,0.01140371],[1.5,3,0.01315938],[2,3,0.0148457],[2.5,3,0.0166198],[3,3,0.01843168],[3.5,3,0.01996274],[4,3,0.02036749],[4.5,3,0.01988819],[5,3,0.01879363],[5.5,3,0.01784356],[6,3,0.01662869],[6.5,3,0.01418338],[7,3,0.01220751],[7.5,3,0.01068294],[8,3,0.00897195],[8.5,3,0.00713635],[9,3,0.005929598],[9.5,3,0.004762196],[10,3,0.003640079],[10.5,3,0.002730215],[11,3,0.001963565],[11.5,3,0.001383902],[12,3,0.000977993],[12.5,3,0.0006551074],[13,3,0.0004047196],[13.5,3,0.0002417829],[14,3,0.0001476967],[14.5,3,9.511718e-005],[15,3,6.420651e-005],[-5,3.5,0.0004548789],[-4.5,3.5,0.0004946627],[-4,3.5,0.0006464385],[-3.5,3.5,0.0009785378],[-3,3.5,0.001444199],[-2.5,3.5,0.002011637],[-2,3.5,0.002678767],[-1.5,3.5,0.003514531],[-1,3.5,0.004598265],[-0.5,3.5,0.005715178],[0,3.5,0.006816225],[0.5,3.5,0.008289365],[1,3.5,0.01013238],[1.5,3.5,0.01178709],[2,3.5,0.01314403],[2.5,3.5,0.01463721],[3,3.5,0.01630229],[3.5,3.5,0.01783916],[4,3.5,0.01821097],[4.5,3.5,0.01773257],[5,3.5,0.01667498],[5.5,3.5,0.01586035],[6,3.5,0.01440971],[6.5,3.5,0.01259961],[7,3.5,0.01113709],[7.5,3.5,0.009736447],[8,3.5,0.008142498],[8.5,3.5,0.006619985],[9,3.5,0.0054469],[9.5,3.5,0.004376777],[10,3.5,0.003387931],[10.5,3.5,0.002567497],[11,3.5,0.001893511],[11.5,3.5,0.001380485],[12,3.5,0.0009531708],[12.5,3.5,0.0006173465],[13,3.5,0.0003837367],[13.5,3.5,0.0002311298],[14,3.5,0.0001398547],[14.5,3.5,8.810128e-005],[15,3.5,5.789307e-005],[-5,4,0.0002346951],[-4.5,4,0.0003390931],[-4,4,0.0005021621],[-3.5,4,0.0007755034],[-3,4,0.00117142],[-2.5,4,0.001638231],[-2,4,0.002179088],[-1.5,4,0.002883999],[-1,4,0.003728086],[-0.5,4,0.004652717],[0,4,0.005667347],[0.5,4,0.00690925],[1,4,0.008371645],[1.5,4,0.009809176],[2,4,0.01108902],[2.5,4,0.01236938],[3,4,0.01366653],[3.5,4,0.01489445],[4,4,0.01541556],[4.5,4,0.01538582],[5,4,0.01433061],[5.5,4,0.01348303],[6,4,0.01243328],[6.5,4,0.01113965],[7,4,0.009905753],[7.5,4,0.008687409],[8,4,0.007212732],[8.5,4,0.00600984],[9,4,0.005055808],[9.5,4,0.003981135],[10,4,0.003111954],[10.5,4,0.002360438],[11,4,0.001745901],[11.5,4,0.001268359],[12,4,0.0008596607],[12.5,4,0.000550915],[13,4,0.0003460691],[13.5,4,0.0002110152],[14,4,0.0001279495],[14.5,4,7.986328e-005],[15,4,5.135033e-005],[-5,4.5,0.0001698713],[-4.5,4.5,0.0002522247],[-4,4.5,0.000373868],[-3.5,4.5,0.0005655993],[-3,4.5,0.000863403],[-2.5,4.5,0.001247074],[-2,4.5,0.001710812],[-1.5,4.5,0.002298051],[-1,4.5,0.00296206],[-0.5,4.5,0.003704701],[0,4.5,0.004559744],[0.5,4.5,0.005666444],[1,4.5,0.006762838],[1.5,4.5,0.007858618],[2,4.5,0.00887587],[2.5,4.5,0.009978087],[3,4.5,0.01095861],[3.5,4.5,0.01189151],[4,4.5,0.0123969],[4.5,4.5,0.01231215],[5,4.5,0.01164148],[5.5,4.5,0.01089155],[6,4.5,0.01022082],[6.5,4.5,0.009588222],[7,4.5,0.008704621],[7.5,4.5,0.007481392],[8,4.5,0.006219859],[8.5,4.5,0.005124464],[9,4.5,0.004226025],[9.5,4.5,0.003408797],[10,4.5,0.002714311],[10.5,4.5,0.002058581],[11,4.5,0.001498668],[11.5,4.5,0.001067243],[12,4.5,0.0007317957],[12.5,4.5,0.0004809838],[13,4.5,0.0003048406],[13.5,4.5,0.0001875519],[14,4.5,0.0001155895],[14.5,4.5,7.332454e-005],[15,4.5,4.726098e-005],[-5,5,0.0001259839],[-4.5,5,0.0001864315],[-4,5,0.00027381],[-3.5,5,0.0004042417],[-3,5,0.000603506],[-2.5,5,0.0008961211],[-2,5,0.001312],[-1.5,5,0.001792428],[-1,5,0.002272279],[-0.5,5,0.002837364],[0,5,0.003548607],[0.5,5,0.00462431],[1,5,0.005253413],[1.5,5,0.006141242],[2,5,0.006924978],[2.5,5,0.007583723],[3,5,0.008362503],[3.5,5,0.009249979],[4,5,0.009770368],[4.5,5,0.009494431],[5,5,0.009127974],[5.5,5,0.008640395],[6,5,0.008217071],[6.5,5,0.007953936],[7,5,0.007102081],[7.5,5,0.006062319],[8,5,0.005172467],[8.5,5,0.004225706],[9,5,0.003536584],[9.5,5,0.002815495],[10,5,0.002233659],[10.5,5,0.001721311],[11,5,0.001262767],[11.5,5,0.0008956209],[12,5,0.0006229865],[12.5,5,0.0004163966],[13,5,0.0002655091],[13.5,5,0.0001655139],[14,5,0.0001050058],[14.5,5,6.889506e-005],[15,5,4.527796e-005],[-5,5.5,9.162851e-005],[-4.5,5.5,0.0001350944],[-4,5.5,0.0001974007],[-3.5,5.5,0.0002883708],[-3,5.5,0.0004234519],[-2.5,5.5,0.0006314362],[-2,5.5,0.0009506683],[-1.5,5.5,0.001308848],[-1,5.5,0.001653377],[-0.5,5.5,0.00208797],[0,5.5,0.00263139],[0.5,5.5,0.003282662],[1,5.5,0.003776289],[1.5,5.5,0.004443085],[2,5.5,0.00508417],[2.5,5.5,0.00559856],[3,5.5,0.006116484],[3.5,5.5,0.006763822],[4,5.5,0.007237205],[4.5,5.5,0.007073841],[5,5.5,0.006930507],[5.5,5.5,0.006719528],[6,5.5,0.006312751],[6.5,5.5,0.005920962],[7,5.5,0.005402633],[7.5,5.5,0.004795589],[8,5.5,0.004102475],[8.5,5.5,0.003352992],[9,5.5,0.002860071],[9.5,5.5,0.002267921],[10,5.5,0.001797238],[10.5,5.5,0.001421605],[11,5.5,0.001061632],[11.5,5.5,0.0007511447],[12,5.5,0.0005267841],[12.5,5.5,0.0003559844],[13,5.5,0.0002294968],[13.5,5.5,0.0001461219],[14,5.5,9.586095e-005],[14.5,5.5,6.488113e-005],[15,5.5,4.317755e-005],[-5,6,6.476804e-005],[-4.5,6,9.521397e-005],[-4,6,0.00013865],[-3.5,6,0.0002013941],[-3,6,0.0002928613],[-2.5,6,0.0004281006],[-2,6,0.0006251941],[-1.5,6,0.0008642646],[-1,6,0.001130087],[-0.5,6,0.001466266],[0,6,0.00187085],[0.5,6,0.002236408],[1,6,0.002587926],[1.5,6,0.003020727],[2,6,0.00349341],[2.5,6,0.003921281],[3,6,0.004296366],[3.5,6,0.004642268],[4,6,0.004909287],[4.5,6,0.004994022],[5,6,0.004933056],[5.5,6,0.00484383],[6,6,0.004563627],[6.5,6,0.004350115],[7,6,0.004100825],[7.5,6,0.003664642],[8,6,0.003080692],[8.5,6,0.002518834],[9,6,0.00209075],[9.5,6,0.001724001],[10,6,0.001382187],[10.5,6,0.00109785],[11,6,0.0008415821],[11.5,6,0.0006155072],[12,6,0.000437278],[12.5,6,0.0002972183],[13,6,0.0001942236],[13.5,6,0.0001269245],[14,6,8.591473e-005],[14.5,6,5.917532e-005],[15,6,3.912643e-005],[-5,6.5,4.425527e-005],[-4.5,6.5,6.488522e-005],[-4,6.5,9.418262e-005],[-3.5,6.5,0.000136066],[-3,6.5,0.0001960708],[-2.5,6.5,0.0002813686],[-2,6.5,0.0003992832],[-1.5,6.5,0.0005517142],[-1,6.5,0.0007396697],[-0.5,6.5,0.0009752643],[0,6.5,0.001255079],[0.5,6.5,0.001512917],[1,6.5,0.001747296],[1.5,6.5,0.002019724],[2,6.5,0.00232418],[2.5,6.5,0.002618556],[3,6.5,0.002883444],[3.5,6.5,0.003133535],[4,6.5,0.003342996],[4.5,6.5,0.003382717],[5,6.5,0.003318907],[5.5,6.5,0.003224797],[6,6.5,0.003156446],[6.5,6.5,0.00309649],[7,6.5,0.002963065],[7.5,6.5,0.002655187],[8,6.5,0.002222098],[8.5,6.5,0.001827538],[9,6.5,0.001520944],[9.5,6.5,0.001253707],[10,6.5,0.001014801],[10.5,6.5,0.0008255297],[11,6.5,0.0006590463],[11.5,6.5,0.0004973785],[12,6.5,0.0003500521],[12.5,6.5,0.0002351663],[13,6.5,0.0001554047],[13.5,6.5,0.0001045233],[14,6.5,7.278701e-005],[14.5,6.5,5.046426e-005],[15,6.5,3.269335e-005],[-5,7,2.911555e-005],[-4.5,7,4.256807e-005],[-4,7,6.156442e-005],[-3.5,7,8.846812e-005],[-3,7,0.0001265047],[-2.5,7,0.000179656],[-2,7,0.0002519964],[-1.5,7,0.0003466832],[-1,7,0.0004666983],[-0.5,7,0.0006181875],[0,7,0.0007996359],[0.5,7,0.0009814304],[1,7,0.001155249],[1.5,7,0.001332416],[2,7,0.001506367],[2.5,7,0.001679478],[3,7,0.00185071],[3.5,7,0.002063192],[4,7,0.002348107],[4.5,7,0.00221149],[5,7,0.002150837],[5.5,7,0.00211658],[6,7,0.002111191],[6.5,7,0.002078967],[7,7,0.001985535],[7.5,7,0.001799088],[8,7,0.001527191],[8.5,7,0.001267901],[9,7,0.001054634],[9.5,7,0.0008692254],[10,7,0.000729815],[10.5,7,0.0006350054],[11,7,0.0005105405],[11.5,7,0.0003726422],[12,7,0.0002545207],[12.5,7,0.0001689197],[13,7,0.0001132909],[13.5,7,7.855338e-005],[14,7,5.603966e-005],[14.5,7,3.887967e-005],[15,7,2.459765e-005],[-5,7.5,1.83765e-005],[-4.5,7.5,2.679906e-005],[-4,7.5,3.861415e-005],[-3.5,7.5,5.522401e-005],[-3,7.5,7.851572e-005],[-2.5,7.5,0.0001107921],[-2,7.5,0.0001544418],[-1.5,7.5,0.0002114922],[-1,7.5,0.0002836878],[-0.5,7.5,0.0003731365],[0,7.5,0.0004796041],[0.5,7.5,0.000599039],[1,7.5,0.0007303341],[1.5,7.5,0.0008498753],[2,7.5,0.0009440681],[2.5,7.5,0.001042782],[3,7.5,0.0011602],[3.5,7.5,0.001296144],[4,7.5,0.001415046],[4.5,7.5,0.001372925],[5,7.5,0.001346269],[5.5,7.5,0.001333717],[6,7.5,0.001333785],[6.5,7.5,0.001316495],[7,7.5,0.001264804],[7.5,7.5,0.001160642],[8,7.5,0.001001858],[8.5,7.5,0.0008396343],[9,7.5,0.0006970846],[9.5,7.5,0.0005722468],[10,7.5,0.0004846482],[10.5,7.5,0.0004259433],[11,7.5,0.0003329455],[11.5,7.5,0.0002349747],[12,7.5,0.0001605743],[12.5,7.5,0.0001082681],[13,7.5,7.421195e-005],[13.5,7.5,5.27031e-005],[14,7.5,3.809068e-005],[14.5,7.5,2.629087e-005],[15,7.5,1.629133e-005],[-5,8,1.109556e-005],[-4.5,8,1.616191e-005],[-4,8,2.322661e-005],[-3.5,8,3.311152e-005],[-3,8,4.692817e-005],[-2.5,8,6.602436e-005],[-2,8,9.177663e-005],[-1.5,8,0.0001252875],[-1,8,0.0001672066],[-0.5,8,0.0002176825],[0,8,0.0002765175],[0.5,8,0.0003458388],[1,8,0.0004268432],[1.5,8,0.0005017466],[2,8,0.0005617295],[2.5,8,0.0006266521],[3,8,0.0007096608],[3.5,8,0.0008003166],[4,8,0.0008464586],[4.5,8,0.0008392808],[5,8,0.0008257861],[5.5,8,0.0008183373],[6,8,0.0008143275],[6.5,8,0.0008085631],[7,8,0.0007913557],[7.5,8,0.0007393145],[8,8,0.0006442568],[8.5,8,0.0005382248],[9,8,0.0004428103],[9.5,8,0.0003584886],[10,8,0.0002902228],[10.5,8,0.0002365976],[11,8,0.0001812033],[11.5,8,0.0001309469],[12,8,9.214342e-005],[12.5,8,6.358873e-005],[13,8,4.413704e-005],[13.5,8,3.146919e-005],[14,8,2.259788e-005],[14.5,8,1.53779e-005],[15,8,9.366106e-006],[-5,8.5,6.399208e-006],[-4.5,8.5,9.334239e-006],[-4,8.5,1.341448e-005],[-3.5,8.5,1.911662e-005],[-3,8.5,2.70925e-005],[-2.5,8.5,3.812565e-005],[-2,8.5,5.298769e-005],[-1.5,8.5,7.222469e-005],[-1,8.5,9.599251e-005],[-0.5,8.5,0.0001240274],[0,8.5,0.0001559371],[0.5,8.5,0.0001923337],[1,8.5,0.0002338707],[1.5,8.5,0.0002763762],[2,8.5,0.0003173889],[2.5,8.5,0.000361462],[3,8.5,0.000412578],[3.5,8.5,0.0004650965],[4,8.5,0.0004946637],[4.5,8.5,0.0004978636],[5,8.5,0.0004952616],[5.5,8.5,0.000492224],[6,8.5,0.0004888952],[6.5,8.5,0.0004887367],[7,8.5,0.0004889782],[7.5,8.5,0.0004662949],[8,8.5,0.0004065263],[8.5,8.5,0.0003346867],[9,8.5,0.0002715891],[9.5,8.5,0.0002168595],[10,8.5,0.0001696181],[10.5,8.5,0.0001306476],[11,8.5,9.834988e-005],[11.5,8.5,7.193484e-005],[12,8.5,5.100255e-005],[12.5,8.5,3.522899e-005],[13,8.5,2.422238e-005],[13.5,8.5,1.690533e-005],[14,8.5,1.179901e-005],[14.5,8.5,7.827571e-006],[15,8.5,4.692644e-006],[-5,9,3.52419e-006],[-4.5,9,5.165248e-006],[-4,9,7.450104e-006],[-3.5,9,1.065146e-005],[-3,9,1.514355e-005],[-2.5,9,2.136756e-005],[-2,9,2.973458e-005],[-1.5,9,4.048974e-005],[-1,9,5.361396e-005],[-0.5,9,6.88455e-005],[0,9,8.586889e-005],[0.5,9,0.0001046809],[1,9,0.0001256308],[1.5,9,0.0001486762],[2,9,0.0001736117],[2.5,9,0.000200589],[3,9,0.0002292169],[3.5,9,0.0002563468],[4,9,0.0002752015],[4.5,9,0.0002841722],[5,9,0.0002875097],[5.5,9,0.000287357],[6,9,0.0002853195],[6.5,9,0.0002847654],[7,9,0.0002855338],[7.5,9,0.0002738103],[8,9,0.0002397165],[8.5,9,0.0001973818],[9,9,0.0001595501],[9.5,9,0.0001263422],[10,9,9.725523e-005],[10.5,9,7.314737e-005],[11,9,5.395888e-005],[11.5,9,3.88833e-005],[12,9,2.719771e-005],[12.5,9,1.84955e-005],[13,9,1.240335e-005],[13.5,9,8.328527e-006],[14,9,5.559846e-006],[14.5,9,3.560558e-006],[15,9,2.099908e-006],[-5,9.5,1.854042e-006],[-4.5,9.5,2.739681e-006],[-4,9.5,3.980201e-006],[-3.5,9.5,5.727809e-006],[-3,9.5,8.189336e-006],[-2.5,9.5,1.160079e-005],[-2,9.5,1.61662e-005],[-1.5,9.5,2.198015e-005],[-1,9.5,2.898201e-005],[-0.5,9.5,3.699161e-005],[0,9.5,4.583762e-005],[0.5,9.5,5.552853e-005],[1,9.5,6.632722e-005],[1.5,9.5,7.858144e-005],[2,9.5,9.245352e-005],[2.5,9.5,0.0001077333],[3,9.5,0.0001235722],[3.5,9.5,0.0001382012],[4,9.5,0.000149505],[4.5,9.5,0.0001566608],[5,9.5,0.0001602828],[5.5,9.5,0.0001611736],[6,9.5,0.0001600916],[6.5,9.5,0.0001581212],[7,9.5,0.0001551826],[7.5,9.5,0.0001468353],[8,9.5,0.0001298971],[8.5,9.5,0.0001090354],[9,9.5,8.883386e-005],[9.5,9.5,7.009047e-005],[10,9.5,5.33809e-005],[10.5,9.5,3.947704e-005],[11,9.5,2.851915e-005],[11.5,9.5,2.011679e-005],[12,9.5,1.37914e-005],[12.5,9.5,9.17878e-006],[13,9.5,5.97337e-006],[13.5,9.5,3.841827e-006],[14,9.5,2.441355e-006],[14.5,9.5,1.502772e-006],[15,9.5,8.709228e-007],[-5,10,9.322845e-007],[-4.5,10,1.392493e-006],[-4,10,2.042776e-006],[-3.5,10,2.964702e-006],[-3,10,4.266762e-006],[-2.5,10,6.067685e-006],[-2,10,8.46114e-006],[-1.5,10,1.147571e-005],[-1,10,1.505818e-005],[-0.5,10,1.910549e-005],[0,10,2.354456e-005],[0.5,10,2.842136e-005],[1,10,3.393007e-005],[1.5,10,4.032199e-005],[2,10,4.771951e-005],[2.5,10,5.594073e-005],[3,10,6.440134e-005],[3.5,10,7.218367e-005],[4,10,7.841497e-005],[4.5,10,8.272415e-005],[5,10,8.524814e-005],[5.5,10,8.624512e-005],[6,10,8.583616e-005],[6.5,10,8.404824e-005],[7,10,8.070927e-005],[7.5,10,7.506846e-005],[8,10,6.676331e-005],[8.5,10,5.687542e-005],[9,10,4.662127e-005],[9.5,10,3.669779e-005],[10,10,2.774117e-005],[10.5,10,2.024683e-005],[11,10,1.435343e-005],[11.5,10,9.904646e-006],[12,10,6.6391e-006],[12.5,10,4.314336e-006],[13,10,2.725526e-006],[13.5,10,1.685681e-006],[14,10,1.024631e-006],[14.5,10,6.076249e-007],[15,10,3.458268e-007],[-5,10.5,4.48198e-007],[-4.5,10.5,6.775336e-007],[-4,10.5,1.004552e-006],[-3.5,10.5,1.470581e-006],[-3,10.5,2.129099e-006],[-2.5,10.5,3.0361e-006],[-2,10.5,4.231632e-006],[-1.5,10.5,5.720954e-006],[-1,10.5,7.47023e-006],[-0.5,10.5,9.42835e-006],[0,10.5,1.157127e-005],[0.5,10.5,1.394542e-005],[1,10.5,1.66741e-005],[1.5,10.5,1.989998e-005],[2,10.5,2.36766e-005],[2.5,10.5,2.787065e-005],[3,10.5,3.214469e-005],[3.5,10.5,3.605619e-005],[4,10.5,3.924458e-005],[4.5,10.5,4.157277e-005],[5,10.5,4.308326e-005],[5.5,10.5,4.382019e-005],[6,10.5,4.371364e-005],[6.5,10.5,4.262382e-005],[7,10.5,4.044648e-005],[7.5,10.5,3.716275e-005],[8,10.5,3.291858e-005],[8.5,10.5,2.806029e-005],[9,10.5,2.29814e-005],[9.5,10.5,1.804319e-005],[10,10.5,1.357717e-005],[10.5,10.5,9.821617e-006],[11,10.5,6.861351e-006],[11.5,10.5,4.643299e-006],[12,10.5,3.044207e-006],[12.5,10.5,1.931266e-006],[13,10.5,1.18694e-006],[13.5,10.5,7.103397e-007],[14,10.5,4.164024e-007],[14.5,10.5,2.391927e-007],[15,10.5,1.336445e-007],[-5,11,2.059511e-007],[-4.5,11,3.150836e-007],[-4,11,4.718769e-007],[-3.5,11,6.959243e-007],[-3,11,1.011857e-006],[-2.5,11,1.444458e-006],[-2,11,2.009965e-006],[-1.5,11,2.707974e-006],[-1,11,3.521077e-006],[-0.5,11,4.427077e-006],[0,11,5.420779e-006],[0.5,11,6.533041e-006],[1,11,7.829953e-006],[1.5,11,9.380772e-006],[2,11,1.120297e-005],[2.5,11,1.321772e-005],[3,11,1.52546e-005],[3.5,11,1.711523e-005],[4,11,1.865752e-005],[4.5,11,1.983423e-005],[5,11,2.065025e-005],[5.5,11,2.108226e-005],[6,11,2.104601e-005],[6.5,11,2.044405e-005],[7,11,1.924188e-005],[7.5,11,1.750343e-005],[8,11,1.537155e-005],[8.5,11,1.302141e-005],[9,11,1.062065e-005],[9.5,11,8.318612e-006],[10,11,6.245282e-006],[10.5,11,4.495787e-006],[11,11,3.110602e-006],[11.5,11,2.074097e-006],[12,11,1.334469e-006],[12.5,11,8.285848e-007],[13,11,4.971585e-007],[13.5,11,2.89598e-007],[14,11,1.648947e-007],[14.5,11,9.219202e-008],[15,11,5.052742e-008],[-5,11.5,9.036716e-008],[-4.5,11.5,1.397819e-007],[-4,11.5,2.111396e-007],[-3.5,11.5,3.130951e-007],[-3,11.5,4.562334e-007],[-2.5,11.5,6.508855e-007],[-2,11.5,9.034493e-007],[-1.5,11.5,1.213215e-006],[-1,11.5,1.572755e-006],[-0.5,11.5,1.973712e-006],[0,11.5,2.416256e-006],[0.5,11.5,2.916534e-006],[1,11.5,3.504874e-006],[1.5,11.5,4.210568e-006],[2,11.5,5.037823e-006],[2.5,11.5,5.948104e-006],[3,11.5,6.865756e-006],[3.5,11.5,7.708427e-006],[4,11.5,8.420829e-006],[4.5,11.5,8.982991e-006],[5,11.5,9.385572e-006],[5.5,11.5,9.598237e-006],[6,11.5,9.566401e-006],[6.5,11.5,9.242875e-006],[7,11.5,8.625235e-006],[7.5,11.5,7.766403e-006],[8,11.5,6.753559e-006],[8.5,11.5,5.676274e-006],[9,11.5,4.607106e-006],[9.5,11.5,3.600966e-006],[10,11.5,2.70153e-006],[10.5,11.5,1.941559e-006],[11,11.5,1.336916e-006],[11.5,11.5,8.83309e-007],[12,11.5,5.608047e-007],[12.5,11.5,3.42521e-007],[13,11.5,2.016732e-007],[13.5,11.5,1.150132e-007],[14,11.5,6.398378e-008],[14.5,11.5,3.494419e-008],[15,11.5,1.876817e-008],[-5,12,3.780181e-008],[-4.5,12,5.903879e-008],[-4,12,8.979193e-008],[-3.5,12,1.336285e-007],[-3,12,1.948034e-007],[-2.5,12,2.773947e-007],[-2,12,3.839061e-007],[-1.5,12,5.14134e-007],[-1,12,6.654761e-007],[-0.5,12,8.35311e-007],[0,12,1.024632e-006],[0.5,12,1.240635e-006],[1,12,1.495516e-006],[1.5,12,1.800203e-006],[2,12,2.155111e-006],[2.5,12,2.543919e-006],[3,12,2.936509e-006],[3.5,12,3.300909e-006],[4,12,3.615147e-006],[4.5,12,3.868309e-006],[5,12,4.049976e-006],[5.5,12,4.139883e-006],[6,12,4.110774e-006],[6.5,12,3.944322e-006],[7,12,3.646451e-006],[7.5,12,3.248903e-006],[8,12,2.796552e-006],[8.5,12,2.33108e-006],[9,12,1.882059e-006],[9.5,12,1.467985e-006],[10,12,1.101535e-006],[10.5,12,7.922234e-007],[11,12,5.450864e-007],[11.5,12,3.587982e-007],[12,12,2.261719e-007],[12.5,12,1.367421e-007],[13,12,7.950654e-008],[13.5,12,4.467075e-008],[14,12,2.44189e-008],[14.5,12,1.307501e-008],[15,12,6.881992e-009],[-5,12.5,1.504492e-008],[-4.5,12.5,2.369478e-008],[-4,12.5,3.623749e-008],[-3.5,12.5,5.405355e-008],[-3,12.5,7.875327e-008],[-2.5,12.5,1.118716e-007],[-2,12.5,1.543908e-007],[-1.5,12.5,2.063738e-007],[-1,12.5,2.67115e-007],[-0.5,12.5,3.360054e-007],[0,12.5,4.1378e-007],[0.5,12.5,5.033094e-007],[1,12.5,6.090056e-007],[1.5,12.5,7.34503e-007],[2,12.5,8.794714e-007],[2.5,12.5,1.037621e-006],[3,12.5,1.197861e-006],[3.5,12.5,1.348323e-006],[4,12.5,1.480048e-006],[4.5,12.5,1.586899e-006],[5,12.5,1.661943e-006],[5.5,12.5,1.694759e-006],[6,12.5,1.673868e-006],[6.5,12.5,1.59345e-006],[7,12.5,1.45885e-006],[7.5,12.5,1.286099e-006],[8,12.5,1.095733e-006],[8.5,12.5,9.055569e-007],[9,12.5,7.269356e-007],[9.5,12.5,5.656225e-007],[10,12.5,4.246139e-007],[10.5,12.5,3.060099e-007],[11,12.5,2.109729e-007],[11.5,12.5,1.389417e-007],[12,12.5,8.742749e-008],[12.5,12.5,5.263994e-008],[13,12.5,3.041392e-008],[13.5,12.5,1.69401e-008],[14,12.5,9.151929e-009],[14.5,12.5,4.825872e-009],[15,12.5,2.493997e-009],[-5,13,5.684774e-009],[-4.5,13,9.021213e-009],[-4,13,1.386449e-008],[-3.5,13,2.072123e-008],[-3,13,3.017103e-008],[-2.5,13,4.276945e-008],[-2,13,5.890053e-008],[-1.5,13,7.866945e-008],[-1,13,1.019623e-007],[-0.5,13,1.287287e-007],[0,13,1.593686e-007],[0.5,13,1.949447e-007],[1,13,2.369293e-007],[1.5,13,2.864132e-007],[2,13,3.430876e-007],[2.5,13,4.046477e-007],[3,13,4.671793e-007],[3.5,13,5.263734e-007],[4,13,5.785747e-007],[4.5,13,6.206837e-007],[5,13,6.491303e-007],[5.5,13,6.593752e-007],[6,13,6.471638e-007],[6.5,13,6.110322e-007],[7,13,5.540929e-007],[7.5,13,4.834934e-007],[8,13,4.077926e-007],[8.5,13,3.340629e-007],[9,13,2.664662e-007],[9.5,13,2.066729e-007],[10,13,1.551423e-007],[10.5,13,1.120652e-007],[11,13,7.75221e-008],[11.5,13,5.121555e-008],[12,13,3.229316e-008],[12.5,13,1.945466e-008],[13,13,1.122776e-008],[13.5,13,6.233282e-009],[14,13,3.346332e-009],[14.5,13,1.74642e-009],[15,13,8.894546e-010],[-5,13.5,2.035356e-009],[-4.5,13.5,3.253851e-009],[-4,13.5,5.026277e-009],[-3.5,13.5,7.53107e-009],[-3,13.5,1.096915e-008],[-2.5,13.5,1.553564e-008],[-2,13.5,2.137815e-008],[-1.5,13.5,2.856811e-008],[-1,13.5,3.712132e-008],[-0.5,13.5,4.708229e-008],[0,13.5,5.863439e-008],[0.5,13.5,7.214904e-008],[1,13.5,8.808913e-008],[1.5,13.5,1.06753e-007],[2,13.5,1.279632e-007],[2.5,13.5,1.50894e-007],[3,13.5,1.741872e-007],[3.5,13.5,1.962971e-007],[4,13.5,2.157783e-007],[4.5,13.5,2.312583e-007],[5,13.5,2.411927e-007],[5.5,13.5,2.43838e-007],[6,13.5,2.377708e-007],[6.5,13.5,2.227418e-007],[7,13.5,2.002118e-007],[7.5,13.5,1.73058e-007],[8,13.5,1.445717e-007],[8.5,13.5,1.173897e-007],[9,13.5,9.297895e-008],[9.5,13.5,7.18053e-008],[10,13.5,5.383756e-008],[10.5,13.5,3.894996e-008],[11,13.5,2.703702e-008],[11.5,13.5,1.793907e-008],[12,13.5,1.13598e-008],[12.5,13.5,6.869026e-009],[13,13.5,3.974916e-009],[13.5,13.5,2.208965e-009],[14,13.5,1.183938e-009],[14.5,13.5,6.14564e-010],[15,13.5,3.099197e-010],[-5,14,6.894894e-010],[-4.5,14,1.110908e-009],[-4,14,1.726374e-009],[-3.5,14,2.596898e-009],[-3,14,3.790426e-009],[-2.5,14,5.374025e-009],[-2,14,7.402798e-009],[-1.5,14,9.912589e-009],[-1,14,1.29262e-008],[-0.5,14,1.647676e-008],[0,14,2.06369e-008],[0.5,14,2.552874e-008],[1,14,3.129023e-008],[1.5,14,3.799551e-008],[2,14,4.55618e-008],[2.5,14,5.369893e-008],[3,14,6.194206e-008],[3.5,14,6.97481e-008],[4,14,7.657658e-008],[4.5,14,8.188927e-008],[5,14,8.510173e-008],[5.5,14,8.560976e-008],[6,14,8.297576e-008],[6.5,14,7.719982e-008],[7,14,6.887156e-008],[7.5,14,5.904885e-008],[8,14,4.890753e-008],[8.5,14,3.937513e-008],[9,14,3.095318e-008],[9.5,14,2.377359e-008],[10,14,1.777666e-008],[10.5,14,1.286342e-008],[11,14,8.952675e-009],[11.5,14,5.965871e-009],[12,14,3.797795e-009],[12.5,14,2.309241e-009],[13,14,1.343307e-009],[13.5,14,7.496267e-010],[14,14,4.026296e-010],[14.5,14,2.087856e-010],[15,14,1.047564e-010],[-5,14.5,2.208027e-010],[-4.5,14.5,3.589015e-010],[-4,14.5,5.61947e-010],[-3.5,14.5,8.503575e-010],[-3,14.5,1.246762e-009],[-2.5,14.5,1.773858e-009],[-2,14.5,2.451479e-009],[-1.5,14.5,3.294666e-009],[-1,14.5,4.315241e-009],[-0.5,14.5,5.527825e-009],[0,14.5,6.957234e-009],[0.5,14.5,8.640315e-009],[1,14.5,1.061534e-008],[1.5,14.5,1.289829e-008],[2,14.5,1.545545e-008],[2.5,14.5,1.818841e-008],[3,14.5,2.094286e-008],[3.5,14.5,2.353565e-008],[4,14.5,2.577876e-008],[4.5,14.5,2.748275e-008],[5,14.5,2.84493e-008],[5.5,14.5,2.848659e-008],[6,14.5,2.746915e-008],[6.5,14.5,2.541871e-008],[7,14.5,2.254606e-008],[7.5,14.5,1.920857e-008],[8,14.5,1.579782e-008],[8.5,14.5,1.262218e-008],[9,14.5,9.848169e-009],[9.5,14.5,7.516127e-009],[10,14.5,5.597004e-009],[10.5,14.5,4.044552e-009],[11,14.5,2.818838e-009],[11.5,14.5,1.885389e-009],[12,14.5,1.206722e-009],[12.5,14.5,7.384929e-010],[13,14.5,4.325175e-010],[13.5,14.5,2.428975e-010],[14,14.5,1.311037e-010],[14.5,14.5,6.815153e-011],[15,14.5,3.416334e-011],[-5,15,6.683681e-011],[-4.5,15,1.097557e-010],[-4,15,1.734784e-010],[-3.5,15,2.647113e-010],[-3,15,3.909017e-010],[-2.5,15,5.596232e-010],[-2,15,7.777459e-010],[-1.5,15,1.050839e-009],[-1,15,1.383444e-009],[-0.5,15,1.780502e-009],[0,15,2.249169e-009],[0.5,15,2.799211e-009],[1,15,3.440019e-009],[1.5,15,4.173971e-009],[2,15,4.988642e-009],[2.5,15,5.852391e-009],[3,15,6.716605e-009],[3.5,15,7.523276e-009],[4,15,8.21217e-009],[4.5,15,8.722865e-009],[5,15,8.994147e-009],[5.5,15,8.9696e-009],[6,15,8.615178e-009],[6.5,15,7.942046e-009],[7,15,7.017922e-009],[7.5,15,5.954162e-009],[8,15,4.872438e-009],[8.5,15,3.869546e-009],[9,15,2.998913e-009],[9.5,15,2.274025e-009],[10,15,1.684882e-009],[10.5,15,1.214291e-009],[11,15,8.463777e-010],[11.5,15,5.676707e-010],[12,15,3.651537e-010],[12.5,15,2.249535e-010],[13,15,1.327439e-010],[13.5,15,7.511369e-011],[14,15,4.081411e-011],[14.5,15,2.131914e-011],[15,15,1.070996e-011]],"colors":[[1,0,0,1]],"centers":[[-5,-5,1.146066e-005],[-4.5,-5,1.79456e-005],[-4,-5,2.826249e-005],[-3.5,-5,4.302606e-005],[-3,-5,6.224708e-005],[-2.5,-5,8.488809e-005],[-2,-5,0.0001095062],[-1.5,-5,0.0001354254],[-1,-5,0.0001627141],[-0.5,-5,0.0001904372],[0,-5,0.0002158209],[0.5,-5,0.0002363462],[1,-5,0.0002515781],[1.5,-5,0.0002620281],[2,-5,0.0002678616],[2.5,-5,0.0002711958],[3,-5,0.0002817815],[3.5,-5,0.0002963529],[4,-5,0.0002775255],[4.5,-5,0.0002341229],[5,-5,0.0002044563],[5.5,-5,0.0001868718],[6,-5,0.0001718011],[6.5,-5,0.0001562788],[7,-5,0.000139487],[7.5,-5,0.0001214974],[8,-5,0.0001024068],[8.5,-5,8.197113e-005],[9,-5,6.098348e-005],[9.5,-5,4.201412e-005],[10,-5,2.748056e-005],[10.5,-5,1.774059e-005],[11,-5,1.156565e-005],[11.5,-5,7.565056e-006],[12,-5,4.852432e-006],[12.5,-5,2.98488e-006],[13,-5,1.735515e-006],[13.5,-5,9.481955e-007],[14,-5,4.875333e-007],[14.5,-5,2.376798e-007],[15,-5,1.111036e-007],[-5,-4.5,2.237209e-005],[-4.5,-4.5,3.180327e-005],[-4,-4.5,4.912843e-005],[-3.5,-4.5,7.398456e-005],[-3,-4.5,0.0001062965],[-2.5,-4.5,0.0001444892],[-2,-4.5,0.0001864878],[-1.5,-4.5,0.0002318203],[-1,-4.5,0.0002811916],[-0.5,-4.5,0.0003324109],[0,-4.5,0.0003786405],[0.5,-4.5,0.0004147704],[1,-4.5,0.0004414355],[1.5,-4.5,0.0004599946],[2,-4.5,0.0004701946],[2.5,-4.5,0.0004742752],[3,-4.5,0.000483124],[3.5,-4.5,0.0004937478],[4,-4.5,0.0004665218],[4.5,-4.5,0.0004111622],[5,-4.5,0.0003687843],[5.5,-4.5,0.0003377087],[6,-4.5,0.0003081178],[6.5,-4.5,0.0002773233],[7,-4.5,0.0002449829],[7.5,-4.5,0.0002116368],[8,-4.5,0.0001777556],[8.5,-4.5,0.0001425613],[9,-4.5,0.0001065677],[9.5,-4.5,7.380409e-005],[10,-4.5,4.854332e-005],[10.5,-4.5,3.157794e-005],[11,-4.5,2.08178e-005],[11.5,-4.5,1.380867e-005],[12,-4.5,8.987108e-006],[12.5,-4.5,5.601876e-006],[13,-4.5,3.293227e-006],[13.5,-4.5,1.814594e-006],[14,-4.5,9.383666e-007],[14.5,-4.5,4.586728e-007],[15,-4.5,2.142398e-007],[-5,-4,4.112635e-005],[-4.5,-4,5.3383e-005],[-4,-4,8.109483e-005],[-3.5,-4,0.0001204543],[-3,-4,0.0001714598],[-2.5,-4,0.0002326409],[-2,-4,0.0003012976],[-1.5,-4,0.000377329],[-1,-4,0.0004629512],[-0.5,-4,0.0005543241],[0,-4,0.0006346446],[0.5,-4,0.0006932709],[1,-4,0.000737445],[1.5,-4,0.0007713854],[2,-4,0.0007912981],[2.5,-4,0.0007962267],[3,-4,0.0007920349],[3.5,-4,0.0007798162],[4,-4,0.0007454995],[4.5,-4,0.0006956958],[5,-4,0.0006462956],[5.5,-4,0.0005930346],[6,-4,0.0005372703],[6.5,-4,0.0004808864],[7,-4,0.0004220904],[7.5,-4,0.0003589195],[8,-4,0.0002956315],[8.5,-4,0.0002345179],[9,-4,0.00017526],[9.5,-4,0.0001223395],[10,-4,8.160929e-005],[10.5,-4,5.40173e-005],[11,-4,3.621635e-005],[11.5,-4,2.437729e-005],[12,-4,1.606446e-005],[12.5,-4,1.012191e-005],[13,-4,6.007861e-006],[13.5,-4,3.339113e-006],[14,-4,1.739856e-006],[14.5,-4,8.55642e-007],[15,-4,4.012541e-007],[-5,-3.5,6.327545e-005],[-4.5,-3.5,8.672065e-005],[-4,-3.5,0.000132443],[-3.5,-3.5,0.0001909186],[-3,-3.5,0.0002642531],[-2.5,-3.5,0.0003558708],[-2,-3.5,0.0004623005],[-1.5,-3.5,0.0005829262],[-1,-3.5,0.000723041],[-0.5,-3.5,0.0008791749],[0,-3.5,0.001015789],[0.5,-3.5,0.001107298],[1,-3.5,0.001176776],[1.5,-3.5,0.001239609],[2,-3.5,0.001280798],[2.5,-3.5,0.001292596],[3,-3.5,0.00128148],[3.5,-3.5,0.001253596],[4,-3.5,0.001213399],[4.5,-3.5,0.001173225],[5,-3.5,0.001118858],[5.5,-3.5,0.001016316],[6,-3.5,0.000907936],[6.5,-3.5,0.0008119827],[7,-3.5,0.0007128508],[7.5,-3.5,0.0005953113],[8,-3.5,0.0004763827],[8.5,-3.5,0.0003724136],[9,-3.5,0.0002785494],[9.5,-3.5,0.0001954714],[10,-3.5,0.0001319688],[10.5,-3.5,8.922184e-005],[11,-3.5,6.112921e-005],[11.5,-3.5,4.182365e-005],[12,-3.5,2.788037e-005],[12.5,-3.5,1.772779e-005],[13,-3.5,1.061292e-005],[13.5,-3.5,5.950761e-006],[14,-3.5,3.129319e-006],[14.5,-3.5,1.553438e-006],[15,-3.5,7.349277e-007],[-5,-3,9.086989e-005],[-4.5,-3,0.0001408084],[-4,-3,0.0002162213],[-3.5,-3,0.0002972678],[-3,-3,0.0003921731],[-2.5,-3,0.0005203511],[-2,-3,0.0006772354],[-1.5,-3,0.0008578757],[-1,-3,0.001066805],[-0.5,-3,0.001303766],[0,-3,0.001525707],[0.5,-3,0.001685078],[1,-3,0.001809782],[1.5,-3,0.001927267],[2,-3,0.002012582],[2.5,-3,0.002050395],[3,-3,0.002030791],[3.5,-3,0.001995545],[4,-3,0.001950376],[4.5,-3,0.00191211],[5,-3,0.001838514],[5.5,-3,0.001657118],[6,-3,0.001466172],[6.5,-3,0.001304361],[7,-3,0.001145016],[7.5,-3,0.0009816453],[8,-3,0.0007547182],[8.5,-3,0.0005843053],[9,-3,0.0004421417],[9.5,-3,0.0003100148],[10,-3,0.0002091989],[10.5,-3,0.0001438895],[11,-3,0.0001008534],[11.5,-3,7.020459e-005],[12,-3,4.734426e-005],[12.5,-3,3.037661e-005],[13,-3,1.834266e-005],[13.5,-3,1.037853e-005],[14,-3,5.51339e-006],[14.5,-3,2.769239e-006],[15,-3,1.327097e-006],[-5,-2.5,0.0001442305],[-4.5,-2.5,0.0002205195],[-4,-2.5,0.0003183437],[-3.5,-2.5,0.0004252557],[-3,-2.5,0.0005536033],[-2.5,-2.5,0.0007333626],[-2,-2.5,0.0009581725],[-1.5,-2.5,0.001216761],[-1,-2.5,0.001502542],[-0.5,-2.5,0.001815174],[0,-2.5,0.002135989],[0.5,-2.5,0.002436238],[1,-2.5,0.002771712],[1.5,-2.5,0.002912784],[2,-2.5,0.003119577],[2.5,-2.5,0.003412267],[3,-2.5,0.003185678],[3.5,-2.5,0.003097689],[4,-2.5,0.003042762],[4.5,-2.5,0.002960254],[5,-2.5,0.002793584],[5.5,-2.5,0.002527546],[6,-2.5,0.002248884],[6.5,-2.5,0.001980725],[7,-2.5,0.001736836],[7.5,-2.5,0.001621885],[8,-2.5,0.001175007],[8.5,-2.5,0.0008900007],[9,-2.5,0.0006791472],[9.5,-2.5,0.0004778789],[10,-2.5,0.0003248813],[10.5,-2.5,0.0002277697],[11,-2.5,0.0001631131],[11.5,-2.5,0.0001154858],[12,-2.5,7.891313e-005],[12.5,-2.5,5.122688e-005],[13,-2.5,3.126673e-005],[13.5,-2.5,1.786575e-005],[14,-2.5,9.590415e-006],[14.5,-2.5,4.880421e-006],[15,-2.5,2.376079e-006],[-5,-2,0.0002165681],[-4.5,-2,0.0003067722],[-4,-2,0.0004096671],[-3.5,-2,0.0005531394],[-3,-2,0.0007481863],[-2.5,-2,0.001009775],[-2,-2,0.001321894],[-1.5,-2,0.001676139],[-1,-2,0.002057783],[-0.5,-2,0.002455546],[0,-2,0.002884016],[0.5,-2,0.003378617],[1,-2,0.00399484],[1.5,-2,0.004209686],[2,-2,0.004477599],[2.5,-2,0.004908534],[3,-2,0.004767854],[3.5,-2,0.004671034],[4,-2,0.004610468],[4.5,-2,0.004455058],[5,-2,0.00412845],[5.5,-2,0.003718831],[6,-2,0.003315366],[6.5,-2,0.002901602],[7,-2,0.002481993],[7.5,-2,0.002087073],[8,-2,0.001655077],[8.5,-2,0.001289009],[9,-2,0.0009762045],[9.5,-2,0.0006985588],[10,-2,0.0004903682],[10.5,-2,0.0003524965],[11,-2,0.000256847],[11.5,-2,0.0001845717],[12,-2,0.0001280628],[12.5,-2,8.446218e-005],[13,-2,5.228769e-005],[13.5,-2,3.0224e-005],[14,-2,1.641317e-005],[14.5,-2,8.478482e-006],[15,-2,4.205202e-006],[-5,-1.5,0.0002572473],[-4.5,-1.5,0.000368599],[-4,-1.5,0.0005108458],[-3.5,-1.5,0.0007113838],[-3,-1.5,0.00098646],[-2.5,-1.5,0.001347442],[-2,-1.5,0.001758092],[-1.5,-1.5,0.002231156],[-1,-1.5,0.002761184],[-0.5,-1.5,0.003279756],[0,-1.5,0.003830162],[0.5,-1.5,0.00450843],[1,-1.5,0.005235973],[1.5,-1.5,0.005755024],[2,-1.5,0.00614239],[2.5,-1.5,0.006588121],[3,-1.5,0.00670407],[3.5,-1.5,0.006766361],[4,-1.5,0.006806706],[4.5,-1.5,0.006441725],[5,-1.5,0.005880353],[5.5,-1.5,0.005301803],[6,-1.5,0.00472252],[6.5,-1.5,0.004121446],[7,-1.5,0.003510275],[7.5,-1.5,0.002946882],[8,-1.5,0.002472274],[8.5,-1.5,0.001842721],[9,-1.5,0.001380014],[9.5,-1.5,0.001006347],[10,-1.5,0.0007270694],[10.5,-1.5,0.0005303404],[11,-1.5,0.000389386],[11.5,-1.5,0.0002825409],[12,-1.5,0.0001988005],[12.5,-1.5,0.0001332903],[13,-1.5,8.380409e-005],[13.5,-1.5,4.912484e-005],[14,-1.5,2.7098e-005],[14.5,-1.5,1.428011e-005],[15,-1.5,7.2475e-006],[-5,-1,0.0002865055],[-4.5,-1,0.0004298274],[-4,-1,0.0006351962],[-3.5,-1,0.0008965567],[-3,-1,0.001237694],[-2.5,-1,0.00169112],[-2,-1,0.002216013],[-1.5,-1,0.002850228],[-1,-1,0.003595983],[-0.5,-1,0.004278233],[0,-1,0.005005321],[0.5,-1,0.005856466],[1,-1,0.006774438],[1.5,-1,0.00756753],[2,-1,0.008159697],[2.5,-1,0.008697816],[3,-1,0.008972279],[3.5,-1,0.009346751],[4,-1,0.009260431],[4.5,-1,0.008641054],[5,-1,0.008031354],[5.5,-1,0.007312033],[6,-1,0.006482272],[6.5,-1,0.005621665],[7,-1,0.004788855],[7.5,-1,0.004025815],[8,-1,0.003372102],[8.5,-1,0.002587227],[9,-1,0.001964825],[9.5,-1,0.001459167],[10,-1,0.001069749],[10.5,-1,0.0007746161],[11,-1,0.0005653928],[11.5,-1,0.0004110206],[12,-1,0.0002910183],[12.5,-1,0.0001969081],[13,-1,0.0001252403],[13.5,-1,7.457702e-005],[14,-1,4.205952e-005],[14.5,-1,2.278823e-005],[15,-1,1.19076e-005],[-5,-0.5,0.0003415506],[-4.5,-0.5,0.0005112388],[-4,-0.5,0.0007535119],[-3.5,-0.5,0.001067786],[-3,-0.5,0.001480715],[-2.5,-0.5,0.00202831],[-2,-0.5,0.002674827],[-1.5,-0.5,0.003443498],[-1,-0.5,0.004380933],[-0.5,-0.5,0.005345846],[0,-0.5,0.006522307],[0.5,-0.5,0.00752635],[1,-0.5,0.008530652],[1.5,-0.5,0.009576255],[2,-0.5,0.01040816],[2.5,-0.5,0.01100683],[3,-0.5,0.01156096],[3.5,-0.5,0.01200045],[4,-0.5,0.01195963],[4.5,-0.5,0.01119982],[5,-0.5,0.0105005],[5.5,-0.5,0.009608634],[6,-0.5,0.0084919],[6.5,-0.5,0.007322086],[7,-0.5,0.006249747],[7.5,-0.5,0.005257108],[8,-0.5,0.004327205],[8.5,-0.5,0.003470191],[9,-0.5,0.002690308],[9.5,-0.5,0.002049894],[10,-0.5,0.001525061],[10.5,-0.5,0.001089799],[11,-0.5,0.0007881323],[11.5,-0.5,0.000570958],[12,-0.5,0.0004018823],[12.5,-0.5,0.0002705456],[13,-0.5,0.0001724871],[13.5,-0.5,0.0001042258],[14,-0.5,6.040592e-005],[14.5,-0.5,3.387763e-005],[15,-0.5,1.833233e-005],[-5,0,0.000404885],[-4.5,0,0.0006008795],[-4,0,0.0008751656],[-3.5,0,0.001247549],[-3,0,0.001736728],[-2.5,0,0.002367277],[-2,0,0.003155175],[-1.5,0,0.004007886],[-1,0,0.005100699],[-0.5,0,0.006445236],[0,0,0.007772099],[0.5,0,0.009231257],[1,0,0.0103401],[1.5,0,0.01150674],[2,0,0.01267288],[2.5,0,0.01361682],[3,0,0.01442646],[3.5,0,0.01490689],[4,0,0.01484609],[4.5,0,0.0139666],[5,0,0.01310446],[5.5,0,0.01197915],[6,0,0.01063787],[6.5,0,0.009245371],[7,0,0.007866147],[7.5,0,0.006665982],[8,0,0.005502604],[8.5,0,0.00439601],[9,0,0.003433482],[9.5,0,0.002644341],[10,0,0.001983844],[10.5,0,0.001441925],[11,0,0.001056572],[11.5,0,0.0007655464],[12,0,0.0005293552],[12.5,0,0.0003487939],[13,0,0.0002205494],[13.5,0,0.0001347026],[14,0,8.026281e-005],[14.5,0,4.670894e-005],[15,0,2.629063e-005],[-5,0.5,0.0004554605],[-4.5,0.5,0.0006792492],[-4,0.5,0.000992506],[-3.5,0.5,0.001417752],[-3,0.5,0.001964358],[-2.5,0.5,0.002657247],[-2,0.5,0.00366049],[-1.5,0.5,0.004595141],[-1,0.5,0.00576894],[-0.5,0.5,0.007459366],[0,0.5,0.00887116],[0.5,0.5,0.01045944],[1,0.5,0.01193677],[1.5,0.5,0.01327685],[2,0.5,0.01474313],[2.5,0.5,0.01608566],[3,0.5,0.0168963],[3.5,0.5,0.01734479],[4,0.5,0.01754291],[4.5,0.5,0.01664309],[5,0.5,0.01560924],[5.5,0.5,0.01441156],[6,0.5,0.0129435],[6.5,0.5,0.0114458],[7,0.5,0.00960853],[7.5,0.5,0.008114891],[8,0.5,0.006699379],[8.5,0.5,0.005356373],[9,0.5,0.004211412],[9.5,0.5,0.003240301],[10,0.5,0.002437766],[10.5,0.5,0.001810906],[11,0.5,0.0013517],[11.5,0.5,0.0009828239],[12,0.5,0.0006644871],[12.5,0.5,0.0004252456],[13,0.5,0.0002654321],[13.5,0.5,0.0001632694],[14,0.5,9.966295e-005],[14.5,0.5,6.011502e-005],[15,0.5,3.534319e-005],[-5,1,0.0004830987],[-4.5,1,0.0007234502],[-4,1,0.001059132],[-3.5,1,0.001505477],[-3,1,0.002079298],[-2.5,1,0.002811247],[-2,1,0.003819107],[-1.5,1,0.00497342],[-1,1,0.006137304],[-0.5,1,0.007744041],[0,1,0.009407871],[0.5,1,0.01132513],[1,1,0.01336012],[1.5,1,0.01477714],[2,1,0.01638014],[2.5,1,0.01806297],[3,1,0.01928366],[3.5,1,0.01978466],[4,1,0.01963364],[4.5,1,0.01880674],[5,1,0.01787322],[5.5,1,0.01652871],[6,1,0.01492306],[6.5,1,0.0134116],[7,1,0.01134544],[7.5,1,0.00951534],[8,1,0.007900327],[8.5,1,0.00638144],[9,1,0.005065174],[9.5,1,0.003906359],[10,1,0.002961425],[10.5,1,0.002188874],[11,1,0.001609593],[11.5,1,0.001159115],[12,1,0.0007749737],[12.5,1,0.0004902227],[13,1,0.0003044707],[13.5,1,0.0001883835],[14,1,0.0001171817],[14.5,1,7.308921e-005],[15,1,4.508184e-005],[-5,1.5,0.0004881235],[-4.5,1.5,0.0007284761],[-4,1.5,0.001065357],[-3.5,1.5,0.001517409],[-3,1.5,0.002097141],[-2.5,1.5,0.002801118],[-2,1.5,0.003797894],[-1.5,1.5,0.004950379],[-1,1.5,0.006141036],[-0.5,1.5,0.007881681],[0,1.5,0.009610836],[0.5,1.5,0.01153261],[1,1.5,0.01378538],[1.5,1.5,0.01543789],[2,1.5,0.01713712],[2.5,1.5,0.01884047],[3,1.5,0.02055656],[3.5,1.5,0.02144197],[4,1.5,0.02136368],[4.5,1.5,0.02052283],[5,1.5,0.01955397],[5.5,1.5,0.01827986],[6,1.5,0.0165303],[6.5,1.5,0.0145778],[7,1.5,0.01255322],[7.5,1.5,0.01068448],[8,1.5,0.00909933],[8.5,1.5,0.007381013],[9,1.5,0.005875655],[9.5,1.5,0.004506598],[10,1.5,0.00346514],[10.5,1.5,0.00251358],[11,1.5,0.001780023],[11.5,1.5,0.001254741],[12,1.5,0.0008444163],[12.5,1.5,0.000541629],[13,1.5,0.0003381055],[13.5,1.5,0.0002097826],[14,1.5,0.0001319739],[14.5,1.5,8.473112e-005],[15,1.5,5.482933e-005],[-5,2,0.0004730257],[-4.5,2,0.000700683],[-4,2,0.001024601],[-3.5,2,0.001507917],[-3,2,0.002236372],[-2.5,2,0.002774193],[-2,2,0.003647106],[-1.5,2,0.004689286],[-1,2,0.00595497],[-0.5,2,0.007754375],[0,2,0.009409109],[0.5,2,0.01114514],[1,2,0.01317165],[1.5,2,0.01518674],[2,2,0.01712975],[2.5,2,0.01882221],[3,2,0.0208555],[3.5,2,0.02175393],[4,2,0.02216771],[4.5,2,0.021446],[5,2,0.02035583],[5.5,2,0.01902198],[6,2,0.01742218],[6.5,2,0.01507541],[7,2,0.01309105],[7.5,2,0.01135501],[8,2,0.009636084],[8.5,2,0.007805859],[9,2,0.006227697],[9.5,2,0.005007782],[10,2,0.003783584],[10.5,2,0.002677933],[11,2,0.001877535],[11.5,2,0.001313652],[12,2,0.0009009013],[12.5,2,0.0005906412],[13,2,0.0003692903],[13.5,2,0.0002274542],[14,2,0.0001430757],[14.5,2,9.355474e-005],[15,2,6.283284e-005],[-5,2.5,0.000436783],[-4.5,2.5,0.0006396043],[-4,2.5,0.0009291347],[-3.5,2.5,0.001357811],[-3,2.5,0.001960082],[-2.5,2.5,0.002600539],[-2,2.5,0.00344066],[-1.5,2.5,0.004441527],[-1,2.5,0.005637107],[-0.5,2.5,0.007220118],[0,2.5,0.008836407],[0.5,2.5,0.01041832],[1,2.5,0.01233387],[1.5,2.5,0.01431841],[2,2.5,0.01620566],[2.5,2.5,0.01812897],[3,2.5,0.01996089],[3.5,2.5,0.02129149],[4,2.5,0.02176238],[4.5,2.5,0.02133545],[5,2.5,0.02015966],[5.5,2.5,0.01897212],[6,2.5,0.0177619],[6.5,2.5,0.01513214],[7,2.5,0.01288958],[7.5,2.5,0.0112167],[8,2.5,0.009383049],[8.5,2.5,0.007482089],[9,2.5,0.006160055],[9.5,2.5,0.005133534],[10,2.5,0.003891281],[10.5,2.5,0.002751673],[11,2.5,0.001946777],[11.5,2.5,0.00136052],[12,2.5,0.0009550565],[12.5,2.5,0.0006391706],[13,2.5,0.0003964631],[13.5,2.5,0.0002396129],[14,2.5,0.0001488367],[14.5,2.5,9.741767e-005],[15,2.5,6.643661e-005],[-5,3,0.0004756195],[-4.5,3,0.0005826598],[-4,3,0.0007925265],[-3.5,3,0.001162311],[-3,3,0.001679428],[-2.5,3,0.002329326],[-2,3,0.003104219],[-1.5,3,0.004049969],[-1,3,0.00525702],[-0.5,3,0.006698369],[0,3,0.007916006],[0.5,3,0.009484243],[1,3,0.01140371],[1.5,3,0.01315938],[2,3,0.0148457],[2.5,3,0.0166198],[3,3,0.01843168],[3.5,3,0.01996274],[4,3,0.02036749],[4.5,3,0.01988819],[5,3,0.01879363],[5.5,3,0.01784356],[6,3,0.01662869],[6.5,3,0.01418338],[7,3,0.01220751],[7.5,3,0.01068294],[8,3,0.00897195],[8.5,3,0.00713635],[9,3,0.005929598],[9.5,3,0.004762196],[10,3,0.003640079],[10.5,3,0.002730215],[11,3,0.001963565],[11.5,3,0.001383902],[12,3,0.000977993],[12.5,3,0.0006551074],[13,3,0.0004047196],[13.5,3,0.0002417829],[14,3,0.0001476967],[14.5,3,9.511718e-005],[15,3,6.420651e-005],[-5,3.5,0.0004548789],[-4.5,3.5,0.0004946627],[-4,3.5,0.0006464385],[-3.5,3.5,0.0009785378],[-3,3.5,0.001444199],[-2.5,3.5,0.002011637],[-2,3.5,0.002678767],[-1.5,3.5,0.003514531],[-1,3.5,0.004598265],[-0.5,3.5,0.005715178],[0,3.5,0.006816225],[0.5,3.5,0.008289365],[1,3.5,0.01013238],[1.5,3.5,0.01178709],[2,3.5,0.01314403],[2.5,3.5,0.01463721],[3,3.5,0.01630229],[3.5,3.5,0.01783916],[4,3.5,0.01821097],[4.5,3.5,0.01773257],[5,3.5,0.01667498],[5.5,3.5,0.01586035],[6,3.5,0.01440971],[6.5,3.5,0.01259961],[7,3.5,0.01113709],[7.5,3.5,0.009736447],[8,3.5,0.008142498],[8.5,3.5,0.006619985],[9,3.5,0.0054469],[9.5,3.5,0.004376777],[10,3.5,0.003387931],[10.5,3.5,0.002567497],[11,3.5,0.001893511],[11.5,3.5,0.001380485],[12,3.5,0.0009531708],[12.5,3.5,0.0006173465],[13,3.5,0.0003837367],[13.5,3.5,0.0002311298],[14,3.5,0.0001398547],[14.5,3.5,8.810128e-005],[15,3.5,5.789307e-005],[-5,4,0.0002346951],[-4.5,4,0.0003390931],[-4,4,0.0005021621],[-3.5,4,0.0007755034],[-3,4,0.00117142],[-2.5,4,0.001638231],[-2,4,0.002179088],[-1.5,4,0.002883999],[-1,4,0.003728086],[-0.5,4,0.004652717],[0,4,0.005667347],[0.5,4,0.00690925],[1,4,0.008371645],[1.5,4,0.009809176],[2,4,0.01108902],[2.5,4,0.01236938],[3,4,0.01366653],[3.5,4,0.01489445],[4,4,0.01541556],[4.5,4,0.01538582],[5,4,0.01433061],[5.5,4,0.01348303],[6,4,0.01243328],[6.5,4,0.01113965],[7,4,0.009905753],[7.5,4,0.008687409],[8,4,0.007212732],[8.5,4,0.00600984],[9,4,0.005055808],[9.5,4,0.003981135],[10,4,0.003111954],[10.5,4,0.002360438],[11,4,0.001745901],[11.5,4,0.001268359],[12,4,0.0008596607],[12.5,4,0.000550915],[13,4,0.0003460691],[13.5,4,0.0002110152],[14,4,0.0001279495],[14.5,4,7.986328e-005],[15,4,5.135033e-005],[-5,4.5,0.0001698713],[-4.5,4.5,0.0002522247],[-4,4.5,0.000373868],[-3.5,4.5,0.0005655993],[-3,4.5,0.000863403],[-2.5,4.5,0.001247074],[-2,4.5,0.001710812],[-1.5,4.5,0.002298051],[-1,4.5,0.00296206],[-0.5,4.5,0.003704701],[0,4.5,0.004559744],[0.5,4.5,0.005666444],[1,4.5,0.006762838],[1.5,4.5,0.007858618],[2,4.5,0.00887587],[2.5,4.5,0.009978087],[3,4.5,0.01095861],[3.5,4.5,0.01189151],[4,4.5,0.0123969],[4.5,4.5,0.01231215],[5,4.5,0.01164148],[5.5,4.5,0.01089155],[6,4.5,0.01022082],[6.5,4.5,0.009588222],[7,4.5,0.008704621],[7.5,4.5,0.007481392],[8,4.5,0.006219859],[8.5,4.5,0.005124464],[9,4.5,0.004226025],[9.5,4.5,0.003408797],[10,4.5,0.002714311],[10.5,4.5,0.002058581],[11,4.5,0.001498668],[11.5,4.5,0.001067243],[12,4.5,0.0007317957],[12.5,4.5,0.0004809838],[13,4.5,0.0003048406],[13.5,4.5,0.0001875519],[14,4.5,0.0001155895],[14.5,4.5,7.332454e-005],[15,4.5,4.726098e-005],[-5,5,0.0001259839],[-4.5,5,0.0001864315],[-4,5,0.00027381],[-3.5,5,0.0004042417],[-3,5,0.000603506],[-2.5,5,0.0008961211],[-2,5,0.001312],[-1.5,5,0.001792428],[-1,5,0.002272279],[-0.5,5,0.002837364],[0,5,0.003548607],[0.5,5,0.00462431],[1,5,0.005253413],[1.5,5,0.006141242],[2,5,0.006924978],[2.5,5,0.007583723],[3,5,0.008362503],[3.5,5,0.009249979],[4,5,0.009770368],[4.5,5,0.009494431],[5,5,0.009127974],[5.5,5,0.008640395],[6,5,0.008217071],[6.5,5,0.007953936],[7,5,0.007102081],[7.5,5,0.006062319],[8,5,0.005172467],[8.5,5,0.004225706],[9,5,0.003536584],[9.5,5,0.002815495],[10,5,0.002233659],[10.5,5,0.001721311],[11,5,0.001262767],[11.5,5,0.0008956209],[12,5,0.0006229865],[12.5,5,0.0004163966],[13,5,0.0002655091],[13.5,5,0.0001655139],[14,5,0.0001050058],[14.5,5,6.889506e-005],[15,5,4.527796e-005],[-5,5.5,9.162851e-005],[-4.5,5.5,0.0001350944],[-4,5.5,0.0001974007],[-3.5,5.5,0.0002883708],[-3,5.5,0.0004234519],[-2.5,5.5,0.0006314362],[-2,5.5,0.0009506683],[-1.5,5.5,0.001308848],[-1,5.5,0.001653377],[-0.5,5.5,0.00208797],[0,5.5,0.00263139],[0.5,5.5,0.003282662],[1,5.5,0.003776289],[1.5,5.5,0.004443085],[2,5.5,0.00508417],[2.5,5.5,0.00559856],[3,5.5,0.006116484],[3.5,5.5,0.006763822],[4,5.5,0.007237205],[4.5,5.5,0.007073841],[5,5.5,0.006930507],[5.5,5.5,0.006719528],[6,5.5,0.006312751],[6.5,5.5,0.005920962],[7,5.5,0.005402633],[7.5,5.5,0.004795589],[8,5.5,0.004102475],[8.5,5.5,0.003352992],[9,5.5,0.002860071],[9.5,5.5,0.002267921],[10,5.5,0.001797238],[10.5,5.5,0.001421605],[11,5.5,0.001061632],[11.5,5.5,0.0007511447],[12,5.5,0.0005267841],[12.5,5.5,0.0003559844],[13,5.5,0.0002294968],[13.5,5.5,0.0001461219],[14,5.5,9.586095e-005],[14.5,5.5,6.488113e-005],[15,5.5,4.317755e-005],[-5,6,6.476804e-005],[-4.5,6,9.521397e-005],[-4,6,0.00013865],[-3.5,6,0.0002013941],[-3,6,0.0002928613],[-2.5,6,0.0004281006],[-2,6,0.0006251941],[-1.5,6,0.0008642646],[-1,6,0.001130087],[-0.5,6,0.001466266],[0,6,0.00187085],[0.5,6,0.002236408],[1,6,0.002587926],[1.5,6,0.003020727],[2,6,0.00349341],[2.5,6,0.003921281],[3,6,0.004296366],[3.5,6,0.004642268],[4,6,0.004909287],[4.5,6,0.004994022],[5,6,0.004933056],[5.5,6,0.00484383],[6,6,0.004563627],[6.5,6,0.004350115],[7,6,0.004100825],[7.5,6,0.003664642],[8,6,0.003080692],[8.5,6,0.002518834],[9,6,0.00209075],[9.5,6,0.001724001],[10,6,0.001382187],[10.5,6,0.00109785],[11,6,0.0008415821],[11.5,6,0.0006155072],[12,6,0.000437278],[12.5,6,0.0002972183],[13,6,0.0001942236],[13.5,6,0.0001269245],[14,6,8.591473e-005],[14.5,6,5.917532e-005],[15,6,3.912643e-005],[-5,6.5,4.425527e-005],[-4.5,6.5,6.488522e-005],[-4,6.5,9.418262e-005],[-3.5,6.5,0.000136066],[-3,6.5,0.0001960708],[-2.5,6.5,0.0002813686],[-2,6.5,0.0003992832],[-1.5,6.5,0.0005517142],[-1,6.5,0.0007396697],[-0.5,6.5,0.0009752643],[0,6.5,0.001255079],[0.5,6.5,0.001512917],[1,6.5,0.001747296],[1.5,6.5,0.002019724],[2,6.5,0.00232418],[2.5,6.5,0.002618556],[3,6.5,0.002883444],[3.5,6.5,0.003133535],[4,6.5,0.003342996],[4.5,6.5,0.003382717],[5,6.5,0.003318907],[5.5,6.5,0.003224797],[6,6.5,0.003156446],[6.5,6.5,0.00309649],[7,6.5,0.002963065],[7.5,6.5,0.002655187],[8,6.5,0.002222098],[8.5,6.5,0.001827538],[9,6.5,0.001520944],[9.5,6.5,0.001253707],[10,6.5,0.001014801],[10.5,6.5,0.0008255297],[11,6.5,0.0006590463],[11.5,6.5,0.0004973785],[12,6.5,0.0003500521],[12.5,6.5,0.0002351663],[13,6.5,0.0001554047],[13.5,6.5,0.0001045233],[14,6.5,7.278701e-005],[14.5,6.5,5.046426e-005],[15,6.5,3.269335e-005],[-5,7,2.911555e-005],[-4.5,7,4.256807e-005],[-4,7,6.156442e-005],[-3.5,7,8.846812e-005],[-3,7,0.0001265047],[-2.5,7,0.000179656],[-2,7,0.0002519964],[-1.5,7,0.0003466832],[-1,7,0.0004666983],[-0.5,7,0.0006181875],[0,7,0.0007996359],[0.5,7,0.0009814304],[1,7,0.001155249],[1.5,7,0.001332416],[2,7,0.001506367],[2.5,7,0.001679478],[3,7,0.00185071],[3.5,7,0.002063192],[4,7,0.002348107],[4.5,7,0.00221149],[5,7,0.002150837],[5.5,7,0.00211658],[6,7,0.002111191],[6.5,7,0.002078967],[7,7,0.001985535],[7.5,7,0.001799088],[8,7,0.001527191],[8.5,7,0.001267901],[9,7,0.001054634],[9.5,7,0.0008692254],[10,7,0.000729815],[10.5,7,0.0006350054],[11,7,0.0005105405],[11.5,7,0.0003726422],[12,7,0.0002545207],[12.5,7,0.0001689197],[13,7,0.0001132909],[13.5,7,7.855338e-005],[14,7,5.603966e-005],[14.5,7,3.887967e-005],[15,7,2.459765e-005],[-5,7.5,1.83765e-005],[-4.5,7.5,2.679906e-005],[-4,7.5,3.861415e-005],[-3.5,7.5,5.522401e-005],[-3,7.5,7.851572e-005],[-2.5,7.5,0.0001107921],[-2,7.5,0.0001544418],[-1.5,7.5,0.0002114922],[-1,7.5,0.0002836878],[-0.5,7.5,0.0003731365],[0,7.5,0.0004796041],[0.5,7.5,0.000599039],[1,7.5,0.0007303341],[1.5,7.5,0.0008498753],[2,7.5,0.0009440681],[2.5,7.5,0.001042782],[3,7.5,0.0011602],[3.5,7.5,0.001296144],[4,7.5,0.001415046],[4.5,7.5,0.001372925],[5,7.5,0.001346269],[5.5,7.5,0.001333717],[6,7.5,0.001333785],[6.5,7.5,0.001316495],[7,7.5,0.001264804],[7.5,7.5,0.001160642],[8,7.5,0.001001858],[8.5,7.5,0.0008396343],[9,7.5,0.0006970846],[9.5,7.5,0.0005722468],[10,7.5,0.0004846482],[10.5,7.5,0.0004259433],[11,7.5,0.0003329455],[11.5,7.5,0.0002349747],[12,7.5,0.0001605743],[12.5,7.5,0.0001082681],[13,7.5,7.421195e-005],[13.5,7.5,5.27031e-005],[14,7.5,3.809068e-005],[14.5,7.5,2.629087e-005],[15,7.5,1.629133e-005],[-5,8,1.109556e-005],[-4.5,8,1.616191e-005],[-4,8,2.322661e-005],[-3.5,8,3.311152e-005],[-3,8,4.692817e-005],[-2.5,8,6.602436e-005],[-2,8,9.177663e-005],[-1.5,8,0.0001252875],[-1,8,0.0001672066],[-0.5,8,0.0002176825],[0,8,0.0002765175],[0.5,8,0.0003458388],[1,8,0.0004268432],[1.5,8,0.0005017466],[2,8,0.0005617295],[2.5,8,0.0006266521],[3,8,0.0007096608],[3.5,8,0.0008003166],[4,8,0.0008464586],[4.5,8,0.0008392808],[5,8,0.0008257861],[5.5,8,0.0008183373],[6,8,0.0008143275],[6.5,8,0.0008085631],[7,8,0.0007913557],[7.5,8,0.0007393145],[8,8,0.0006442568],[8.5,8,0.0005382248],[9,8,0.0004428103],[9.5,8,0.0003584886],[10,8,0.0002902228],[10.5,8,0.0002365976],[11,8,0.0001812033],[11.5,8,0.0001309469],[12,8,9.214342e-005],[12.5,8,6.358873e-005],[13,8,4.413704e-005],[13.5,8,3.146919e-005],[14,8,2.259788e-005],[14.5,8,1.53779e-005],[15,8,9.366106e-006],[-5,8.5,6.399208e-006],[-4.5,8.5,9.334239e-006],[-4,8.5,1.341448e-005],[-3.5,8.5,1.911662e-005],[-3,8.5,2.70925e-005],[-2.5,8.5,3.812565e-005],[-2,8.5,5.298769e-005],[-1.5,8.5,7.222469e-005],[-1,8.5,9.599251e-005],[-0.5,8.5,0.0001240274],[0,8.5,0.0001559371],[0.5,8.5,0.0001923337],[1,8.5,0.0002338707],[1.5,8.5,0.0002763762],[2,8.5,0.0003173889],[2.5,8.5,0.000361462],[3,8.5,0.000412578],[3.5,8.5,0.0004650965],[4,8.5,0.0004946637],[4.5,8.5,0.0004978636],[5,8.5,0.0004952616],[5.5,8.5,0.000492224],[6,8.5,0.0004888952],[6.5,8.5,0.0004887367],[7,8.5,0.0004889782],[7.5,8.5,0.0004662949],[8,8.5,0.0004065263],[8.5,8.5,0.0003346867],[9,8.5,0.0002715891],[9.5,8.5,0.0002168595],[10,8.5,0.0001696181],[10.5,8.5,0.0001306476],[11,8.5,9.834988e-005],[11.5,8.5,7.193484e-005],[12,8.5,5.100255e-005],[12.5,8.5,3.522899e-005],[13,8.5,2.422238e-005],[13.5,8.5,1.690533e-005],[14,8.5,1.179901e-005],[14.5,8.5,7.827571e-006],[15,8.5,4.692644e-006],[-5,9,3.52419e-006],[-4.5,9,5.165248e-006],[-4,9,7.450104e-006],[-3.5,9,1.065146e-005],[-3,9,1.514355e-005],[-2.5,9,2.136756e-005],[-2,9,2.973458e-005],[-1.5,9,4.048974e-005],[-1,9,5.361396e-005],[-0.5,9,6.88455e-005],[0,9,8.586889e-005],[0.5,9,0.0001046809],[1,9,0.0001256308],[1.5,9,0.0001486762],[2,9,0.0001736117],[2.5,9,0.000200589],[3,9,0.0002292169],[3.5,9,0.0002563468],[4,9,0.0002752015],[4.5,9,0.0002841722],[5,9,0.0002875097],[5.5,9,0.000287357],[6,9,0.0002853195],[6.5,9,0.0002847654],[7,9,0.0002855338],[7.5,9,0.0002738103],[8,9,0.0002397165],[8.5,9,0.0001973818],[9,9,0.0001595501],[9.5,9,0.0001263422],[10,9,9.725523e-005],[10.5,9,7.314737e-005],[11,9,5.395888e-005],[11.5,9,3.88833e-005],[12,9,2.719771e-005],[12.5,9,1.84955e-005],[13,9,1.240335e-005],[13.5,9,8.328527e-006],[14,9,5.559846e-006],[14.5,9,3.560558e-006],[15,9,2.099908e-006],[-5,9.5,1.854042e-006],[-4.5,9.5,2.739681e-006],[-4,9.5,3.980201e-006],[-3.5,9.5,5.727809e-006],[-3,9.5,8.189336e-006],[-2.5,9.5,1.160079e-005],[-2,9.5,1.61662e-005],[-1.5,9.5,2.198015e-005],[-1,9.5,2.898201e-005],[-0.5,9.5,3.699161e-005],[0,9.5,4.583762e-005],[0.5,9.5,5.552853e-005],[1,9.5,6.632722e-005],[1.5,9.5,7.858144e-005],[2,9.5,9.245352e-005],[2.5,9.5,0.0001077333],[3,9.5,0.0001235722],[3.5,9.5,0.0001382012],[4,9.5,0.000149505],[4.5,9.5,0.0001566608],[5,9.5,0.0001602828],[5.5,9.5,0.0001611736],[6,9.5,0.0001600916],[6.5,9.5,0.0001581212],[7,9.5,0.0001551826],[7.5,9.5,0.0001468353],[8,9.5,0.0001298971],[8.5,9.5,0.0001090354],[9,9.5,8.883386e-005],[9.5,9.5,7.009047e-005],[10,9.5,5.33809e-005],[10.5,9.5,3.947704e-005],[11,9.5,2.851915e-005],[11.5,9.5,2.011679e-005],[12,9.5,1.37914e-005],[12.5,9.5,9.17878e-006],[13,9.5,5.97337e-006],[13.5,9.5,3.841827e-006],[14,9.5,2.441355e-006],[14.5,9.5,1.502772e-006],[15,9.5,8.709228e-007],[-5,10,9.322845e-007],[-4.5,10,1.392493e-006],[-4,10,2.042776e-006],[-3.5,10,2.964702e-006],[-3,10,4.266762e-006],[-2.5,10,6.067685e-006],[-2,10,8.46114e-006],[-1.5,10,1.147571e-005],[-1,10,1.505818e-005],[-0.5,10,1.910549e-005],[0,10,2.354456e-005],[0.5,10,2.842136e-005],[1,10,3.393007e-005],[1.5,10,4.032199e-005],[2,10,4.771951e-005],[2.5,10,5.594073e-005],[3,10,6.440134e-005],[3.5,10,7.218367e-005],[4,10,7.841497e-005],[4.5,10,8.272415e-005],[5,10,8.524814e-005],[5.5,10,8.624512e-005],[6,10,8.583616e-005],[6.5,10,8.404824e-005],[7,10,8.070927e-005],[7.5,10,7.506846e-005],[8,10,6.676331e-005],[8.5,10,5.687542e-005],[9,10,4.662127e-005],[9.5,10,3.669779e-005],[10,10,2.774117e-005],[10.5,10,2.024683e-005],[11,10,1.435343e-005],[11.5,10,9.904646e-006],[12,10,6.6391e-006],[12.5,10,4.314336e-006],[13,10,2.725526e-006],[13.5,10,1.685681e-006],[14,10,1.024631e-006],[14.5,10,6.076249e-007],[15,10,3.458268e-007],[-5,10.5,4.48198e-007],[-4.5,10.5,6.775336e-007],[-4,10.5,1.004552e-006],[-3.5,10.5,1.470581e-006],[-3,10.5,2.129099e-006],[-2.5,10.5,3.0361e-006],[-2,10.5,4.231632e-006],[-1.5,10.5,5.720954e-006],[-1,10.5,7.47023e-006],[-0.5,10.5,9.42835e-006],[0,10.5,1.157127e-005],[0.5,10.5,1.394542e-005],[1,10.5,1.66741e-005],[1.5,10.5,1.989998e-005],[2,10.5,2.36766e-005],[2.5,10.5,2.787065e-005],[3,10.5,3.214469e-005],[3.5,10.5,3.605619e-005],[4,10.5,3.924458e-005],[4.5,10.5,4.157277e-005],[5,10.5,4.308326e-005],[5.5,10.5,4.382019e-005],[6,10.5,4.371364e-005],[6.5,10.5,4.262382e-005],[7,10.5,4.044648e-005],[7.5,10.5,3.716275e-005],[8,10.5,3.291858e-005],[8.5,10.5,2.806029e-005],[9,10.5,2.29814e-005],[9.5,10.5,1.804319e-005],[10,10.5,1.357717e-005],[10.5,10.5,9.821617e-006],[11,10.5,6.861351e-006],[11.5,10.5,4.643299e-006],[12,10.5,3.044207e-006],[12.5,10.5,1.931266e-006],[13,10.5,1.18694e-006],[13.5,10.5,7.103397e-007],[14,10.5,4.164024e-007],[14.5,10.5,2.391927e-007],[15,10.5,1.336445e-007],[-5,11,2.059511e-007],[-4.5,11,3.150836e-007],[-4,11,4.718769e-007],[-3.5,11,6.959243e-007],[-3,11,1.011857e-006],[-2.5,11,1.444458e-006],[-2,11,2.009965e-006],[-1.5,11,2.707974e-006],[-1,11,3.521077e-006],[-0.5,11,4.427077e-006],[0,11,5.420779e-006],[0.5,11,6.533041e-006],[1,11,7.829953e-006],[1.5,11,9.380772e-006],[2,11,1.120297e-005],[2.5,11,1.321772e-005],[3,11,1.52546e-005],[3.5,11,1.711523e-005],[4,11,1.865752e-005],[4.5,11,1.983423e-005],[5,11,2.065025e-005],[5.5,11,2.108226e-005],[6,11,2.104601e-005],[6.5,11,2.044405e-005],[7,11,1.924188e-005],[7.5,11,1.750343e-005],[8,11,1.537155e-005],[8.5,11,1.302141e-005],[9,11,1.062065e-005],[9.5,11,8.318612e-006],[10,11,6.245282e-006],[10.5,11,4.495787e-006],[11,11,3.110602e-006],[11.5,11,2.074097e-006],[12,11,1.334469e-006],[12.5,11,8.285848e-007],[13,11,4.971585e-007],[13.5,11,2.89598e-007],[14,11,1.648947e-007],[14.5,11,9.219202e-008],[15,11,5.052742e-008],[-5,11.5,9.036716e-008],[-4.5,11.5,1.397819e-007],[-4,11.5,2.111396e-007],[-3.5,11.5,3.130951e-007],[-3,11.5,4.562334e-007],[-2.5,11.5,6.508855e-007],[-2,11.5,9.034493e-007],[-1.5,11.5,1.213215e-006],[-1,11.5,1.572755e-006],[-0.5,11.5,1.973712e-006],[0,11.5,2.416256e-006],[0.5,11.5,2.916534e-006],[1,11.5,3.504874e-006],[1.5,11.5,4.210568e-006],[2,11.5,5.037823e-006],[2.5,11.5,5.948104e-006],[3,11.5,6.865756e-006],[3.5,11.5,7.708427e-006],[4,11.5,8.420829e-006],[4.5,11.5,8.982991e-006],[5,11.5,9.385572e-006],[5.5,11.5,9.598237e-006],[6,11.5,9.566401e-006],[6.5,11.5,9.242875e-006],[7,11.5,8.625235e-006],[7.5,11.5,7.766403e-006],[8,11.5,6.753559e-006],[8.5,11.5,5.676274e-006],[9,11.5,4.607106e-006],[9.5,11.5,3.600966e-006],[10,11.5,2.70153e-006],[10.5,11.5,1.941559e-006],[11,11.5,1.336916e-006],[11.5,11.5,8.83309e-007],[12,11.5,5.608047e-007],[12.5,11.5,3.42521e-007],[13,11.5,2.016732e-007],[13.5,11.5,1.150132e-007],[14,11.5,6.398378e-008],[14.5,11.5,3.494419e-008],[15,11.5,1.876817e-008],[-5,12,3.780181e-008],[-4.5,12,5.903879e-008],[-4,12,8.979193e-008],[-3.5,12,1.336285e-007],[-3,12,1.948034e-007],[-2.5,12,2.773947e-007],[-2,12,3.839061e-007],[-1.5,12,5.14134e-007],[-1,12,6.654761e-007],[-0.5,12,8.35311e-007],[0,12,1.024632e-006],[0.5,12,1.240635e-006],[1,12,1.495516e-006],[1.5,12,1.800203e-006],[2,12,2.155111e-006],[2.5,12,2.543919e-006],[3,12,2.936509e-006],[3.5,12,3.300909e-006],[4,12,3.615147e-006],[4.5,12,3.868309e-006],[5,12,4.049976e-006],[5.5,12,4.139883e-006],[6,12,4.110774e-006],[6.5,12,3.944322e-006],[7,12,3.646451e-006],[7.5,12,3.248903e-006],[8,12,2.796552e-006],[8.5,12,2.33108e-006],[9,12,1.882059e-006],[9.5,12,1.467985e-006],[10,12,1.101535e-006],[10.5,12,7.922234e-007],[11,12,5.450864e-007],[11.5,12,3.587982e-007],[12,12,2.261719e-007],[12.5,12,1.367421e-007],[13,12,7.950654e-008],[13.5,12,4.467075e-008],[14,12,2.44189e-008],[14.5,12,1.307501e-008],[15,12,6.881992e-009],[-5,12.5,1.504492e-008],[-4.5,12.5,2.369478e-008],[-4,12.5,3.623749e-008],[-3.5,12.5,5.405355e-008],[-3,12.5,7.875327e-008],[-2.5,12.5,1.118716e-007],[-2,12.5,1.543908e-007],[-1.5,12.5,2.063738e-007],[-1,12.5,2.67115e-007],[-0.5,12.5,3.360054e-007],[0,12.5,4.1378e-007],[0.5,12.5,5.033094e-007],[1,12.5,6.090056e-007],[1.5,12.5,7.34503e-007],[2,12.5,8.794714e-007],[2.5,12.5,1.037621e-006],[3,12.5,1.197861e-006],[3.5,12.5,1.348323e-006],[4,12.5,1.480048e-006],[4.5,12.5,1.586899e-006],[5,12.5,1.661943e-006],[5.5,12.5,1.694759e-006],[6,12.5,1.673868e-006],[6.5,12.5,1.59345e-006],[7,12.5,1.45885e-006],[7.5,12.5,1.286099e-006],[8,12.5,1.095733e-006],[8.5,12.5,9.055569e-007],[9,12.5,7.269356e-007],[9.5,12.5,5.656225e-007],[10,12.5,4.246139e-007],[10.5,12.5,3.060099e-007],[11,12.5,2.109729e-007],[11.5,12.5,1.389417e-007],[12,12.5,8.742749e-008],[12.5,12.5,5.263994e-008],[13,12.5,3.041392e-008],[13.5,12.5,1.69401e-008],[14,12.5,9.151929e-009],[14.5,12.5,4.825872e-009],[15,12.5,2.493997e-009],[-5,13,5.684774e-009],[-4.5,13,9.021213e-009],[-4,13,1.386449e-008],[-3.5,13,2.072123e-008],[-3,13,3.017103e-008],[-2.5,13,4.276945e-008],[-2,13,5.890053e-008],[-1.5,13,7.866945e-008],[-1,13,1.019623e-007],[-0.5,13,1.287287e-007],[0,13,1.593686e-007],[0.5,13,1.949447e-007],[1,13,2.369293e-007],[1.5,13,2.864132e-007],[2,13,3.430876e-007],[2.5,13,4.046477e-007],[3,13,4.671793e-007],[3.5,13,5.263734e-007],[4,13,5.785747e-007],[4.5,13,6.206837e-007],[5,13,6.491303e-007],[5.5,13,6.593752e-007],[6,13,6.471638e-007],[6.5,13,6.110322e-007],[7,13,5.540929e-007],[7.5,13,4.834934e-007],[8,13,4.077926e-007],[8.5,13,3.340629e-007],[9,13,2.664662e-007],[9.5,13,2.066729e-007],[10,13,1.551423e-007],[10.5,13,1.120652e-007],[11,13,7.75221e-008],[11.5,13,5.121555e-008],[12,13,3.229316e-008],[12.5,13,1.945466e-008],[13,13,1.122776e-008],[13.5,13,6.233282e-009],[14,13,3.346332e-009],[14.5,13,1.74642e-009],[15,13,8.894546e-010],[-5,13.5,2.035356e-009],[-4.5,13.5,3.253851e-009],[-4,13.5,5.026277e-009],[-3.5,13.5,7.53107e-009],[-3,13.5,1.096915e-008],[-2.5,13.5,1.553564e-008],[-2,13.5,2.137815e-008],[-1.5,13.5,2.856811e-008],[-1,13.5,3.712132e-008],[-0.5,13.5,4.708229e-008],[0,13.5,5.863439e-008],[0.5,13.5,7.214904e-008],[1,13.5,8.808913e-008],[1.5,13.5,1.06753e-007],[2,13.5,1.279632e-007],[2.5,13.5,1.50894e-007],[3,13.5,1.741872e-007],[3.5,13.5,1.962971e-007],[4,13.5,2.157783e-007],[4.5,13.5,2.312583e-007],[5,13.5,2.411927e-007],[5.5,13.5,2.43838e-007],[6,13.5,2.377708e-007],[6.5,13.5,2.227418e-007],[7,13.5,2.002118e-007],[7.5,13.5,1.73058e-007],[8,13.5,1.445717e-007],[8.5,13.5,1.173897e-007],[9,13.5,9.297895e-008],[9.5,13.5,7.18053e-008],[10,13.5,5.383756e-008],[10.5,13.5,3.894996e-008],[11,13.5,2.703702e-008],[11.5,13.5,1.793907e-008],[12,13.5,1.13598e-008],[12.5,13.5,6.869026e-009],[13,13.5,3.974916e-009],[13.5,13.5,2.208965e-009],[14,13.5,1.183938e-009],[14.5,13.5,6.14564e-010],[15,13.5,3.099197e-010],[-5,14,6.894894e-010],[-4.5,14,1.110908e-009],[-4,14,1.726374e-009],[-3.5,14,2.596898e-009],[-3,14,3.790426e-009],[-2.5,14,5.374025e-009],[-2,14,7.402798e-009],[-1.5,14,9.912589e-009],[-1,14,1.29262e-008],[-0.5,14,1.647676e-008],[0,14,2.06369e-008],[0.5,14,2.552874e-008],[1,14,3.129023e-008],[1.5,14,3.799551e-008],[2,14,4.55618e-008],[2.5,14,5.369893e-008],[3,14,6.194206e-008],[3.5,14,6.97481e-008],[4,14,7.657658e-008],[4.5,14,8.188927e-008],[5,14,8.510173e-008],[5.5,14,8.560976e-008],[6,14,8.297576e-008],[6.5,14,7.719982e-008],[7,14,6.887156e-008],[7.5,14,5.904885e-008],[8,14,4.890753e-008],[8.5,14,3.937513e-008],[9,14,3.095318e-008],[9.5,14,2.377359e-008],[10,14,1.777666e-008],[10.5,14,1.286342e-008],[11,14,8.952675e-009],[11.5,14,5.965871e-009],[12,14,3.797795e-009],[12.5,14,2.309241e-009],[13,14,1.343307e-009],[13.5,14,7.496267e-010],[14,14,4.026296e-010],[14.5,14,2.087856e-010],[15,14,1.047564e-010],[-5,14.5,2.208027e-010],[-4.5,14.5,3.589015e-010],[-4,14.5,5.61947e-010],[-3.5,14.5,8.503575e-010],[-3,14.5,1.246762e-009],[-2.5,14.5,1.773858e-009],[-2,14.5,2.451479e-009],[-1.5,14.5,3.294666e-009],[-1,14.5,4.315241e-009],[-0.5,14.5,5.527825e-009],[0,14.5,6.957234e-009],[0.5,14.5,8.640315e-009],[1,14.5,1.061534e-008],[1.5,14.5,1.289829e-008],[2,14.5,1.545545e-008],[2.5,14.5,1.818841e-008],[3,14.5,2.094286e-008],[3.5,14.5,2.353565e-008],[4,14.5,2.577876e-008],[4.5,14.5,2.748275e-008],[5,14.5,2.84493e-008],[5.5,14.5,2.848659e-008],[6,14.5,2.746915e-008],[6.5,14.5,2.541871e-008],[7,14.5,2.254606e-008],[7.5,14.5,1.920857e-008],[8,14.5,1.579782e-008],[8.5,14.5,1.262218e-008],[9,14.5,9.848169e-009],[9.5,14.5,7.516127e-009],[10,14.5,5.597004e-009],[10.5,14.5,4.044552e-009],[11,14.5,2.818838e-009],[11.5,14.5,1.885389e-009],[12,14.5,1.206722e-009],[12.5,14.5,7.384929e-010],[13,14.5,4.325175e-010],[13.5,14.5,2.428975e-010],[14,14.5,1.311037e-010],[14.5,14.5,6.815153e-011],[15,14.5,3.416334e-011],[-5,15,6.683681e-011],[-4.5,15,1.097557e-010],[-4,15,1.734784e-010],[-3.5,15,2.647113e-010],[-3,15,3.909017e-010],[-2.5,15,5.596232e-010],[-2,15,7.777459e-010],[-1.5,15,1.050839e-009],[-1,15,1.383444e-009],[-0.5,15,1.780502e-009],[0,15,2.249169e-009],[0.5,15,2.799211e-009],[1,15,3.440019e-009],[1.5,15,4.173971e-009],[2,15,4.988642e-009],[2.5,15,5.852391e-009],[3,15,6.716605e-009],[3.5,15,7.523276e-009],[4,15,8.21217e-009],[4.5,15,8.722865e-009],[5,15,8.994147e-009],[5.5,15,8.9696e-009],[6,15,8.615178e-009],[6.5,15,7.942046e-009],[7,15,7.017922e-009],[7.5,15,5.954162e-009],[8,15,4.872438e-009],[8.5,15,3.869546e-009],[9,15,2.998913e-009],[9.5,15,2.274025e-009],[10,15,1.684882e-009],[10.5,15,1.214291e-009],[11,15,8.463777e-010],[11.5,15,5.676707e-010],[12,15,3.651537e-010],[12.5,15,2.249535e-010],[13,15,1.327439e-010],[13.5,15,7.511369e-011],[14,15,4.081411e-011],[14.5,15,2.131914e-011],[15,15,1.070996e-011]],"ignoreExtent":false,"flags":128},"175":{"id":175,"type":"linestrip","material":{"lit":false},"vertices":[[-5,-5,1.022358e-005],[-4.5,-5,1.103743e-005],[-4,-5,1.318789e-005],[-3.5,-5,1.733092e-005],[-3,-5,2.381306e-005],[-2.5,-5,3.314449e-005],[-2,-5,4.359418e-005],[-1.5,-5,5.19521e-005],[-1,-5,6.044171e-005],[-0.5,-5,7.313126e-005],[0,-5,8.640253e-005],[0.5,-5,0.0001005138],[1,-5,0.0001102444],[1.5,-5,0.0001154121],[2,-5,0.0001190983],[2.5,-5,0.0001233736],[3,-5,0.0001256147],[3.5,-5,0.00012259],[4,-5,0.0001157894],[4.5,-5,0.0001086386],[5,-5,0.0001049353],[5.5,-5,8.828189e-005],[6,-5,7.543493e-005],[6.5,-5,6.864825e-005],[7,-5,5.840666e-005],[7.5,-5,4.45579e-005],[8,-5,3.183114e-005],[8.5,-5,2.234906e-005],[9,-5,1.559137e-005],[9.5,-5,1.020451e-005],[10,-5,5.836081e-006],[10.5,-5,2.841148e-006],[11,-5,1.225685e-006],[11.5,-5,6.253013e-007],[12,-5,5.327451e-007],[12.5,-5,3.598875e-007],[13,-5,1.156378e-007],[13.5,-5,2.102437e-008],[14,-5,5.226523e-009],[14.5,-5,3.040046e-009],[15,-5,2.85077e-009],[-5,-4.5,2.069899e-005],[-4.5,-4.5,2.323213e-005],[-4,-4.5,2.796962e-005],[-3.5,-4.5,3.652863e-005],[-3,-4.5,4.949964e-005],[-2.5,-4.5,6.759587e-005],[-2,-4.5,8.810722e-005],[-1.5,-4.5,0.0001060953],[-1,-4.5,0.0001253808],[-0.5,-4.5,0.0001508964],[0,-4.5,0.0001732802],[0.5,-4.5,0.0001988866],[1,-4.5,0.0002193458],[1.5,-4.5,0.0002361014],[2,-4.5,0.0002513785],[2.5,-4.5,0.000254854],[3,-4.5,0.0002701528],[3.5,-4.5,0.0002624239],[4,-4.5,0.0002550682],[4.5,-4.5,0.000227719],[5,-4.5,0.0002172898],[5.5,-4.5,0.0001897597],[6,-4.5,0.0001542671],[6.5,-4.5,0.0001306565],[7,-4.5,0.0001076691],[7.5,-4.5,8.379408e-005],[8,-4.5,6.203192e-005],[8.5,-4.5,4.501045e-005],[9,-4.5,3.251125e-005],[9.5,-4.5,2.241555e-005],[10,-4.5,1.39816e-005],[10.5,-4.5,7.73611e-006],[11,-4.5,3.969961e-006],[11.5,-4.5,2.757731e-006],[12,-4.5,3.067434e-006],[12.5,-4.5,2.218886e-006],[13,-4.5,7.071512e-007],[13.5,-4.5,1.201952e-007],[14,-4.5,2.868395e-008],[14.5,-4.5,1.734266e-008],[15,-4.5,1.631938e-008],[-5,-4,3.352584e-005],[-4.5,-4,4.049542e-005],[-4,-4,5.214016e-005],[-3.5,-4,7.164631e-005],[-3,-4,9.75771e-005],[-2.5,-4,0.0001290569],[-2,-4,0.0001648147],[-1.5,-4,0.0002035393],[-1,-4,0.0002423821],[-0.5,-4,0.000277658],[0,-4,0.0003232002],[0.5,-4,0.0003687174],[1,-4,0.0004159762],[1.5,-4,0.0004651498],[2,-4,0.0004928176],[2.5,-4,0.0004928212],[3,-4,0.0005251947],[3.5,-4,0.0005273412],[4,-4,0.0005226753],[4.5,-4,0.0004599563],[5,-4,0.0004064105],[5.5,-4,0.0003559696],[6,-4,0.0003009211],[6.5,-4,0.0002474041],[7,-4,0.0002011052],[7.5,-4,0.0001592811],[8,-4,0.000121623],[8.5,-4,9.098976e-005],[9,-4,6.744826e-005],[9.5,-4,4.823787e-005],[10,-4,3.217364e-005],[10.5,-4,1.960548e-005],[11,-4,1.114271e-005],[11.5,-4,7.964444e-006],[12,-4,8.288611e-006],[12.5,-4,5.854569e-006],[13,-4,1.943029e-006],[13.5,-4,3.967965e-007],[14,-4,1.23337e-007],[14.5,-4,7.886126e-008],[15,-4,7.274446e-008],[-5,-3.5,4.694144e-005],[-4.5,-3.5,6.258528e-005],[-4,-3.5,9.101545e-005],[-3.5,-3.5,0.00013642],[-3,-3.5,0.0001868602],[-2.5,-3.5,0.0002384381],[-2,-3.5,0.0002940826],[-1.5,-3.5,0.0003653084],[-1,-3.5,0.0004399769],[-0.5,-3.5,0.0005080014],[0,-3.5,0.0005885305],[0.5,-3.5,0.0006682887],[1,-3.5,0.0007498874],[1.5,-3.5,0.000833279],[2,-3.5,0.0008924222],[2.5,-3.5,0.0009123718],[3,-3.5,0.0009428606],[3.5,-3.5,0.0009573893],[4,-3.5,0.0009415169],[4.5,-3.5,0.0008548081],[5,-3.5,0.0007581883],[5.5,-3.5,0.0006717056],[6,-3.5,0.0005796183],[6.5,-3.5,0.0004669309],[7,-3.5,0.0003775536],[7.5,-3.5,0.0003008256],[8,-3.5,0.0002325377],[8.5,-3.5,0.0001762972],[9,-3.5,0.0001313529],[9.5,-3.5,9.535612e-005],[10,-3.5,6.698941e-005],[10.5,-3.5,4.461935e-005],[11,-3.5,2.761806e-005],[11.5,-3.5,1.714831e-005],[12,-3.5,1.286095e-005],[12.5,-3.5,8.092175e-006],[13,-3.5,3.111339e-006],[13.5,-3.5,9.944935e-007],[14,-3.5,4.363978e-007],[14.5,-3.5,2.895069e-007],[15,-3.5,2.559674e-007],[-5,-3,6.283291e-005],[-4.5,-3,9.251746e-005],[-4,-3,0.0001471142],[-3.5,-3,0.0002297969],[-3,-3,0.0003121307],[-2.5,-3,0.0004001023],[-2,-3,0.0005044413],[-1.5,-3,0.0006330485],[-1,-3,0.0007576831],[-0.5,-3,0.0008751507],[0,-3,0.001006732],[0.5,-3,0.001166384],[1,-3,0.001326532],[1.5,-3,0.001486583],[2,-3,0.001629656],[2.5,-3,0.001650302],[3,-3,0.001670886],[3.5,-3,0.001664893],[4,-3,0.001608517],[4.5,-3,0.001487143],[5,-3,0.001343528],[5.5,-3,0.001183867],[6,-3,0.001019062],[6.5,-3,0.0008525469],[7,-3,0.0006903602],[7.5,-3,0.0005486452],[8,-3,0.0004263885],[8.5,-3,0.0003392127],[9,-3,0.0002492083],[9.5,-3,0.0001753216],[10,-3,0.0001275737],[10.5,-3,9.32021e-005],[11,-3,6.687154e-005],[11.5,-3,3.787612e-005],[12,-3,2.088358e-005],[12.5,-3,1.170257e-005],[13,-3,5.556998e-006],[13.5,-3,2.539993e-006],[14,-3,1.321974e-006],[14.5,-3,8.805005e-007],[15,-3,7.269303e-007],[-5,-2.5,8.317419e-005],[-4.5,-2.5,0.0001287157],[-4,-2.5,0.0002019633],[-3.5,-2.5,0.0003048374],[-3,-2.5,0.0004271723],[-2.5,-2.5,0.0005840205],[-2,-2.5,0.0007828551],[-1.5,-2.5,0.001003632],[-1,-2.5,0.0012005],[-0.5,-2.5,0.001396805],[0,-2.5,0.001652539],[0.5,-2.5,0.001935821],[1,-2.5,0.002161358],[1.5,-2.5,0.002413663],[2,-2.5,0.002598478],[2.5,-2.5,0.002700389],[3,-2.5,0.002729656],[3.5,-2.5,0.002745303],[4,-2.5,0.002672856],[4.5,-2.5,0.002461533],[5,-2.5,0.002237034],[5.5,-2.5,0.001982393],[6,-2.5,0.001717946],[6.5,-2.5,0.001445418],[7,-2.5,0.001193795],[7.5,-2.5,0.0009578428],[8,-2.5,0.0007434991],[8.5,-2.5,0.0005971643],[9,-2.5,0.0004381583],[9.5,-2.5,0.0003055494],[10,-2.5,0.0002236774],[10.5,-2.5,0.0001677915],[11,-2.5,0.0001164451],[11.5,-2.5,6.919775e-005],[12,-2.5,4.006087e-005],[12.5,-2.5,2.291202e-005],[13,-2.5,1.204298e-005],[13.5,-2.5,6.205087e-006],[14,-2.5,3.464526e-006],[14.5,-2.5,2.271179e-006],[15,-2.5,1.70555e-006],[-5,-2,0.0001093828],[-4.5,-2,0.0001718271],[-4,-2,0.0002713594],[-3.5,-2,0.0004102671],[-3,-2,0.000606637],[-2.5,-2,0.0008770715],[-2,-2,0.001139808],[-1.5,-2,0.001483015],[-1,-2,0.001788957],[-0.5,-2,0.002087827],[0,-2,0.002536748],[0.5,-2,0.0029928],[1,-2,0.003341079],[1.5,-2,0.003675518],[2,-2,0.003913896],[2.5,-2,0.004107998],[3,-2,0.004225537],[3.5,-2,0.004248867],[4,-2,0.004191169],[4.5,-2,0.003934869],[5,-2,0.003611415],[5.5,-2,0.003176386],[6,-2,0.002791955],[6.5,-2,0.002387549],[7,-2,0.001954116],[7.5,-2,0.001578444],[8,-2,0.001234379],[8.5,-2,0.0009568912],[9,-2,0.0007223485],[9.5,-2,0.0005203434],[10,-2,0.0003706378],[10.5,-2,0.0002685387],[11,-2,0.0001847422],[11.5,-2,0.0001216623],[12,-2,7.370643e-005],[12.5,-2,4.366912e-005],[13,-2,2.452365e-005],[13.5,-2,1.356168e-005],[14,-2,7.896796e-006],[14.5,-2,4.996446e-006],[15,-2,3.35533e-006],[-5,-1.5,0.0001523374],[-4.5,-1.5,0.0002392302],[-4,-1.5,0.0003742393],[-3.5,-1.5,0.0005713367],[-3,-1.5,0.0008743363],[-2.5,-1.5,0.001272012],[-2,-1.5,0.001629646],[-1.5,-1.5,0.002085713],[-1,-1.5,0.002547823],[-0.5,-1.5,0.003040961],[0,-1.5,0.003670539],[0.5,-1.5,0.004358381],[1,-1.5,0.004943821],[1.5,-1.5,0.00544283],[2,-1.5,0.00580303],[2.5,-1.5,0.006083431],[3,-1.5,0.006255389],[3.5,-1.5,0.006263445],[4,-1.5,0.006231996],[4.5,-1.5,0.005965018],[5,-1.5,0.005470508],[5.5,-1.5,0.004870108],[6,-1.5,0.004405541],[6.5,-1.5,0.003792814],[7,-1.5,0.003072205],[7.5,-1.5,0.002482503],[8,-1.5,0.001972175],[8.5,-1.5,0.001516801],[9,-1.5,0.001173589],[9.5,-1.5,0.0008787551],[10,-1.5,0.0006026304],[10.5,-1.5,0.0004110658],[11,-1.5,0.0002887338],[11.5,-1.5,0.0002069897],[12,-1.5,0.0001298236],[12.5,-1.5,7.481874e-005],[13,-1.5,4.448937e-005],[13.5,-1.5,2.616374e-005],[14,-1.5,1.557472e-005],[14.5,-1.5,9.368237e-006],[15,-1.5,5.62355e-006],[-5,-1,0.0002136888],[-4.5,-1,0.0003333431],[-4,-1,0.0005013258],[-3.5,-1,0.0007607086],[-3,-1,0.001157601],[-2.5,-1,0.001646307],[-2,-1,0.002114639],[-1.5,-1,0.002712109],[-1,-1,0.003365271],[-0.5,-1,0.004141956],[0,-1,0.00507151],[0.5,-1,0.006020527],[1,-1,0.006887507],[1.5,-1,0.007648244],[2,-1,0.008223148],[2.5,-1,0.008626256],[3,-1,0.008855464],[3.5,-1,0.008871207],[4,-1,0.00884754],[4.5,-1,0.008500398],[5,-1,0.007864795],[5.5,-1,0.007073918],[6,-1,0.006330473],[6.5,-1,0.005483005],[7,-1,0.004553836],[7.5,-1,0.003709831],[8,-1,0.002964412],[8.5,-1,0.002274409],[9,-1,0.001699988],[9.5,-1,0.001265323],[10,-1,0.0008921861],[10.5,-1,0.0006048291],[11,-1,0.0004253611],[11.5,-1,0.000295162],[12,-1,0.0001942911],[12.5,-1,0.0001176085],[13,-1,7.219736e-005],[13.5,-1,4.442909e-005],[14,-1,2.67895e-005],[14.5,-1,1.538571e-005],[15,-1,8.441748e-006],[-5,-0.5,0.0002770707],[-4.5,-0.5,0.0004275369],[-4,-0.5,0.000636605],[-3.5,-0.5,0.0009603131],[-3,-0.5,0.001425674],[-2.5,-0.5,0.001933799],[-2,-0.5,0.002548859],[-1.5,-0.5,0.003309109],[-1,-0.5,0.004211568],[-0.5,-0.5,0.005276743],[0,-0.5,0.006549423],[0.5,-0.5,0.007879692],[1,-0.5,0.009108986],[1.5,-0.5,0.01014582],[2,-0.5,0.01098782],[2.5,-0.5,0.01164978],[3,-0.5,0.01205101],[3.5,-0.5,0.01212336],[4,-0.5,0.01200153],[4.5,-0.5,0.01161393],[5,-0.5,0.0109102],[5.5,-0.5,0.009806607],[6,-0.5,0.008608765],[6.5,-0.5,0.007440581],[7,-0.5,0.006348599],[7.5,-0.5,0.00516997],[8,-0.5,0.004149278],[8.5,-0.5,0.003169118],[9,-0.5,0.002349548],[9.5,-0.5,0.001728037],[10,-0.5,0.001260719],[10.5,-0.5,0.0008905967],[11,-0.5,0.0006411481],[11.5,-0.5,0.0004184243],[12,-0.5,0.000272737],[12.5,-0.5,0.0001730736],[13,-0.5,0.0001066119],[13.5,-0.5,6.724011e-005],[14,-0.5,4.109306e-005],[14.5,-0.5,2.316816e-005],[15,-0.5,1.219354e-005],[-5,0,0.0003209224],[-4.5,0,0.0005134235],[-4,0,0.0007811526],[-3.5,0,0.001155516],[-3,0,0.001656024],[-2.5,0,0.002254551],[-2,0,0.003031094],[-1.5,0,0.003888661],[-1,0,0.005040784],[-0.5,0,0.006334898],[0,0,0.007765371],[0.5,0,0.009466299],[1,0,0.01119199],[1.5,0,0.01256562],[2,0,0.01372529],[2.5,0,0.01473876],[3,0,0.01540706],[3.5,0,0.01564884],[4,0,0.01544846],[4.5,0,0.01491693],[5,0,0.01399857],[5.5,0,0.01274496],[6,0,0.0112606],[6.5,0,0.009763407],[7,0,0.008276966],[7.5,0,0.006787405],[8,0,0.005504527],[8.5,0,0.004269429],[9,0,0.003221611],[9.5,0,0.002376397],[10,0,0.001716942],[10.5,0,0.00123333],[11,0,0.0008741664],[11.5,0,0.0005769234],[12,0,0.0003682042],[12.5,0,0.0002346456],[13,0,0.0001459708],[13.5,0,9.256292e-005],[14,0,5.688789e-005],[14.5,0,3.211615e-005],[15,0,1.696517e-005],[-5,0.5,0.0003593504],[-4.5,0.5,0.0005729354],[-4,0.5,0.0008808575],[-3.5,0.5,0.001294177],[-3,0.5,0.001845379],[-2.5,0.5,0.002540605],[-2,0.5,0.003430831],[-1.5,0.5,0.004458918],[-1,0.5,0.00577658],[-0.5,0.5,0.007278748],[0,0.5,0.00888849],[0.5,0.5,0.0108683],[1,0.5,0.01292731],[1.5,0.5,0.01464605],[2,0.5,0.01622435],[2.5,0.5,0.0176887],[3,0.5,0.01859451],[3.5,0.5,0.01897665],[4,0.5,0.01891237],[4.5,0.5,0.01822779],[5,0.5,0.01715785],[5.5,0.5,0.01579985],[6,0.5,0.01400261],[6.5,0.5,0.01221664],[7,0.5,0.01039986],[7.5,0.5,0.008624972],[8,0.5,0.006918869],[8.5,0.5,0.00535449],[9,0.5,0.004114065],[9.5,0.5,0.003094596],[10,0.5,0.002246772],[10.5,0.5,0.001618595],[11,0.5,0.001131231],[11.5,0.5,0.0007361925],[12,0.5,0.000475811],[12.5,0.5,0.0003090613],[13,0.5,0.0001923105],[13.5,0.5,0.0001255382],[14,0.5,7.590792e-005],[14.5,0.5,4.042645e-005],[15,0.5,2.15269e-005],[-5,1,0.0003643982],[-4.5,1,0.0005795374],[-4,1,0.0008874646],[-3.5,1,0.001318627],[-3,1,0.001903539],[-2.5,1,0.002640061],[-2,1,0.003648189],[-1.5,1,0.004888502],[-1,1,0.006412107],[-0.5,1,0.00805085],[0,1,0.009901511],[0.5,1,0.01209675],[1,1,0.01435261],[1.5,1,0.01647073],[2,1,0.01847414],[2.5,1,0.0200619],[3,1,0.0211264],[3.5,1,0.0216377],[4,1,0.0216244],[4.5,1,0.02106634],[5,1,0.01994348],[5.5,1,0.0184497],[6,1,0.0165868],[6.5,1,0.01455675],[7,1,0.01253109],[7.5,1,0.01031779],[8,1,0.008187769],[8.5,1,0.006428245],[9,1,0.005023892],[9.5,1,0.003807913],[10,1,0.002770599],[10.5,1,0.002007507],[11,1,0.00140066],[11.5,1,0.0009044391],[12,1,0.0005990935],[12.5,1,0.0003827609],[13,1,0.0002416837],[13.5,1,0.0001558393],[14,1,9.121562e-005],[14.5,1,4.713896e-005],[15,1,2.500576e-005],[-5,1.5,0.0003455108],[-4.5,1.5,0.0005551279],[-4,1.5,0.0008628383],[-3.5,1.5,0.001296526],[-3,1.5,0.00186034],[-2.5,1.5,0.002663968],[-2,1.5,0.00371794],[-1.5,1.5,0.005032463],[-1,1.5,0.006645836],[-0.5,1.5,0.008366022],[0,1.5,0.01037911],[0.5,1.5,0.01270083],[1,1.5,0.01515518],[1.5,1.5,0.01752976],[2,1.5,0.01964824],[2.5,1.5,0.02122781],[3,1.5,0.02255295],[3.5,1.5,0.02335948],[4,1.5,0.02360849],[4.5,1.5,0.02303916],[5,1.5,0.02184367],[5.5,1.5,0.02030874],[6,1.5,0.01852154],[6.5,1.5,0.01637102],[7,1.5,0.01390242],[7.5,1.5,0.01151755],[8,1.5,0.009351893],[8.5,1.5,0.007391411],[9,1.5,0.005814216],[9.5,1.5,0.004422875],[10,1.5,0.003235861],[10.5,1.5,0.002316762],[11,1.5,0.001612895],[11.5,1.5,0.001057941],[12,1.5,0.0007024868],[12.5,1.5,0.0004524757],[13,1.5,0.0002974331],[13.5,1.5,0.0001774191],[14,1.5,9.890705e-005],[14.5,1.5,5.325396e-005],[15,1.5,2.803842e-005],[-5,2,0.0003409614],[-4.5,2,0.0005307131],[-4,2,0.0008160742],[-3.5,2,0.001263079],[-3,2,0.0018135],[-2.5,2,0.002597569],[-2,2,0.003628985],[-1.5,2,0.004940581],[-1,2,0.006500918],[-0.5,2,0.008184075],[0,2,0.01021786],[0.5,2,0.01263407],[1,2,0.01503569],[1.5,2,0.01736504],[2,2,0.01951774],[2.5,2,0.0212557],[3,2,0.02268321],[3.5,2,0.02363331],[4,2,0.02414834],[4.5,2,0.02376101],[5,2,0.02266241],[5.5,2,0.0212651],[6,2,0.01941289],[6.5,2,0.01724573],[7,2,0.01475658],[7.5,2,0.01238161],[8,2,0.01026321],[8.5,2,0.008163632],[9,2,0.006462851],[9.5,2,0.004850509],[10,2,0.003579124],[10.5,2,0.002564469],[11,2,0.001775951],[11.5,2,0.001188788],[12,2,0.0007933318],[12.5,2,0.0005130534],[13,2,0.0003411908],[13.5,2,0.0002005631],[14,2,0.0001106596],[14.5,2,5.930436e-005],[15,2,3.103482e-005],[-5,2.5,0.0003000559],[-4.5,2.5,0.0004768276],[-4,2.5,0.0007216992],[-3.5,2.5,0.001124901],[-3,2.5,0.001645392],[-2.5,2.5,0.002346262],[-2,2.5,0.003320437],[-1.5,2.5,0.004548607],[-1,2.5,0.006013018],[-0.5,2.5,0.007680389],[0,2.5,0.009529932],[0.5,2.5,0.01184772],[1,2.5,0.01408372],[1.5,2.5,0.01629614],[2,2.5,0.01836063],[2.5,2.5,0.02033126],[3,2.5,0.02178055],[3.5,2.5,0.02278999],[4,2.5,0.02336491],[4.5,2.5,0.02317152],[5,2.5,0.02239033],[5.5,2.5,0.02120073],[6,2.5,0.01938601],[6.5,2.5,0.01719356],[7,2.5,0.01495233],[7.5,2.5,0.01264051],[8,2.5,0.01051631],[8.5,2.5,0.008530149],[9,2.5,0.006660428],[9.5,2.5,0.005018294],[10,2.5,0.003767823],[10.5,2.5,0.002747339],[11,2.5,0.001898165],[11.5,2.5,0.001264519],[12,2.5,0.0008248388],[12.5,2.5,0.0005535765],[13,2.5,0.0003488683],[13.5,2.5,0.0002137795],[14,2.5,0.0001238781],[14.5,2.5,6.509617e-005],[15,2.5,3.39542e-005],[-5,3,0.0002451456],[-4.5,3,0.0003979953],[-4,3,0.0006155684],[-3.5,3,0.0009329339],[-3,3,0.001345857],[-2.5,3,0.001988185],[-2,3,0.002849144],[-1.5,3,0.003884822],[-1,3,0.005174107],[-0.5,3,0.006660471],[0,3,0.008364709],[0.5,3,0.01031511],[1,3,0.01235795],[1.5,3,0.01450223],[2,3,0.0165354],[2.5,3,0.01844088],[3,3,0.01983436],[3.5,3,0.02082883],[4,3,0.02133358],[4.5,3,0.02128823],[5,3,0.02070817],[5.5,3,0.01979863],[6,3,0.01821609],[6.5,3,0.01626259],[7,3,0.01424764],[7.5,3,0.01213086],[8,3,0.01019035],[8.5,3,0.008268859],[9,3,0.006466097],[9.5,3,0.004957418],[10,3,0.003724719],[10.5,3,0.002745169],[11,3,0.001940411],[11.5,3,0.001308244],[12,3,0.0008651852],[12.5,3,0.0005758551],[13,3,0.0003675276],[13.5,3,0.0002385012],[14,3,0.0001492127],[14.5,3,7.314443e-005],[15,3,3.72384e-005],[-5,3.5,0.000191473],[-4.5,3.5,0.0003055675],[-4,3.5,0.0004799721],[-3.5,3.5,0.0007385293],[-3,3.5,0.00105315],[-2.5,3.5,0.001554785],[-2,3.5,0.002266804],[-1.5,3.5,0.003137288],[-1,3.5,0.004213018],[-0.5,3.5,0.00551402],[0,3.5,0.006971536],[0.5,3.5,0.008633199],[1,3.5,0.01031888],[1.5,3.5,0.01222394],[2,3.5,0.01411774],[2.5,3.5,0.01564666],[3,3.5,0.01703367],[3.5,3.5,0.0180793],[4,3.5,0.01851347],[4.5,3.5,0.01859006],[5,3.5,0.01817553],[5.5,3.5,0.01741717],[6,3.5,0.01626178],[6.5,3.5,0.01468736],[7,3.5,0.01292671],[7.5,3.5,0.01101026],[8,3.5,0.00922345],[8.5,3.5,0.00756762],[9,3.5,0.005978955],[9.5,3.5,0.004642918],[10,3.5,0.003535577],[10.5,3.5,0.002570006],[11,3.5,0.001824409],[11.5,3.5,0.001259681],[12,3.5,0.0008608194],[12.5,3.5,0.0005828783],[13,3.5,0.0003787368],[13.5,3.5,0.0002434085],[14,3.5,0.0001505519],[14.5,3.5,7.745292e-005],[15,3.5,4.095705e-005],[-5,4,0.0001357248],[-4.5,4,0.0002255649],[-4,4,0.0003600576],[-3.5,4,0.0005506782],[-3,4,0.0008093445],[-2.5,4,0.001201118],[-2,4,0.001722471],[-1.5,4,0.002406857],[-1,4,0.003227973],[-0.5,4,0.004253789],[0,4,0.005446336],[0.5,4,0.00679932],[1,4,0.008197852],[1.5,4,0.009746222],[2,4,0.01123769],[2.5,4,0.01263365],[3,4,0.01386307],[3.5,4,0.01486972],[4,4,0.01537166],[4.5,4,0.01555655],[5,4,0.01531142],[5.5,4,0.01466002],[6,4,0.01369118],[6.5,4,0.01242437],[7,4,0.01104347],[7.5,4,0.00943411],[8,4,0.007913624],[8.5,4,0.006556359],[9,4,0.005241971],[9.5,4,0.004096225],[10,4,0.00312608],[10.5,4,0.002266139],[11,4,0.001634431],[11.5,4,0.001181684],[12,4,0.0008147332],[12.5,4,0.0005360313],[13,4,0.0003551107],[13.5,4,0.0002228951],[14,4,0.0001345286],[14.5,4,8.036088e-005],[15,4,4.527891e-005],[-5,4.5,9.175488e-005],[-4.5,4.5,0.0001702452],[-4,4.5,0.0002787684],[-3.5,4.5,0.0004011895],[-3,4.5,0.0006036633],[-2.5,4.5,0.0009055028],[-2,4.5,0.001253258],[-1.5,4.5,0.001756835],[-1,4.5,0.002347103],[-0.5,4.5,0.003127987],[0,4.5,0.004042109],[0.5,4.5,0.005106617],[1,4.5,0.006192939],[1.5,4.5,0.007418597],[2,4.5,0.008585056],[2.5,4.5,0.009734198],[3,4.5,0.01068319],[3.5,4.5,0.01160795],[4,4.5,0.01216311],[4.5,4.5,0.01230431],[5,4.5,0.01219629],[5.5,4.5,0.0116391],[6,4.5,0.01086461],[6.5,4.5,0.009935535],[7,4.5,0.008893002],[7.5,4.5,0.007675618],[8,4.5,0.006469627],[8.5,4.5,0.005335872],[9,4.5,0.004239196],[9.5,4.5,0.003336417],[10,4.5,0.002581238],[10.5,4.5,0.00194901],[11,4.5,0.001429508],[11.5,4.5,0.001017616],[12,4.5,0.0006895923],[12.5,4.5,0.0004562329],[13,4.5,0.0003049306],[13.5,4.5,0.0001927053],[14,4.5,0.0001205748],[14.5,4.5,8.014286e-005],[15,4.5,4.763472e-005],[-5,5,6.032499e-005],[-4.5,5,0.000114748],[-4,5,0.0001958683],[-3.5,5,0.0002798341],[-3,5,0.0004227639],[-2.5,5,0.0006159834],[-2,5,0.0008340368],[-1.5,5,0.001175063],[-1,5,0.001604252],[-0.5,5,0.002135911],[0,5,0.002823031],[0.5,5,0.003610313],[1,5,0.004428428],[1.5,5,0.005360837],[2,5,0.006242996],[2.5,5,0.007136788],[3,5,0.007872204],[3.5,5,0.008579962],[4,5,0.009138712],[4.5,5,0.009257336],[5,5,0.009088407],[5.5,5,0.008715331],[6,5,0.008201539],[6.5,5,0.007474039],[7,5,0.006758425],[7.5,5,0.005975006],[8,5,0.005082452],[8.5,5,0.00420611],[9,5,0.003310479],[9.5,5,0.00262351],[10,5,0.002034788],[10.5,5,0.001534017],[11,5,0.001145809],[11.5,5,0.0008252413],[12,5,0.0005840532],[12.5,5,0.0003841919],[13,5,0.00025616],[13.5,5,0.0001605711],[14,5,0.0001009728],[14.5,5,6.846455e-005],[15,5,4.156762e-005],[-5,5.5,4.071459e-005],[-4.5,5.5,7.248987e-005],[-4,5.5,0.0001175548],[-3.5,5.5,0.0001754333],[-3,5.5,0.0002725037],[-2.5,5.5,0.0003915074],[-2,5.5,0.0005407106],[-1.5,5.5,0.0007630697],[-1,5.5,0.001073102],[-0.5,5.5,0.001467081],[0,5.5,0.001899732],[0.5,5.5,0.002422495],[1,5.5,0.003023389],[1.5,5.5,0.003678903],[2,5.5,0.004317995],[2.5,5.5,0.004929435],[3,5.5,0.005439891],[3.5,5.5,0.005974792],[4,5.5,0.006449226],[4.5,5.5,0.006568059],[5,5.5,0.006392343],[5.5,5.5,0.006178277],[6,5.5,0.005841099],[6.5,5.5,0.005366211],[7,5.5,0.00491371],[7.5,5.5,0.00440813],[8,5.5,0.003769309],[8.5,5.5,0.003128091],[9,5.5,0.002493047],[9.5,5.5,0.001977746],[10,5.5,0.001534838],[10.5,5.5,0.001155969],[11,5.5,0.0008337336],[11.5,5.5,0.0006096122],[12,5.5,0.0004446622],[12.5,5.5,0.000298758],[13,5.5,0.0001974886],[13.5,5.5,0.0001262647],[14,5.5,8.007965e-005],[14.5,5.5,5.097986e-005],[15,5.5,2.836026e-005],[-5,6,2.813023e-005],[-4.5,6,4.718639e-005],[-4,6,7.148542e-005],[-3.5,6,0.0001038899],[-3,6,0.0001559349],[-2.5,6,0.0002291988],[-2,6,0.0003405188],[-1.5,6,0.0004787663],[-1,6,0.0006712947],[-0.5,6,0.0009144391],[0,6,0.001201808],[0.5,6,0.001546483],[1,6,0.001954929],[1.5,6,0.002375743],[2,6,0.002777802],[2.5,6,0.003198453],[3,6,0.003560571],[3.5,6,0.003904834],[4,6,0.004191277],[4.5,6,0.004318669],[5,6,0.004246304],[5.5,6,0.004124766],[6,6,0.003921967],[6.5,6,0.003677133],[7,6,0.003448613],[7.5,6,0.003106221],[8,6,0.002672201],[8.5,6,0.002201206],[9,6,0.001761596],[9.5,6,0.001407846],[10,6,0.001043773],[10.5,6,0.0008115911],[11,6,0.0006150946],[11.5,6,0.0004439898],[12,6,0.0003138883],[12.5,6,0.0002098406],[13,6,0.0001371921],[13.5,6,9.075213e-005],[14,6,5.848986e-005],[14.5,6,3.420993e-005],[15,6,1.692253e-005],[-5,6.5,1.825595e-005],[-4.5,6.5,3.017732e-005],[-4,6.5,4.384929e-005],[-3.5,6.5,6.186057e-005],[-3,6.5,9.250914e-005],[-2.5,6.5,0.0001328797],[-2,6.5,0.0001933983],[-1.5,6.5,0.0002760926],[-1,6.5,0.0003815318],[-0.5,6.5,0.000523039],[0,6.5,0.0007330539],[0.5,6.5,0.0009513506],[1,6.5,0.001200232],[1.5,6.5,0.001447563],[2,6.5,0.001689532],[2.5,6.5,0.001970515],[3,6.5,0.002186827],[3.5,6.5,0.002399571],[4,6.5,0.002558368],[4.5,6.5,0.002626491],[5,6.5,0.002610845],[5.5,6.5,0.002588239],[6,6.5,0.00248109],[6.5,6.5,0.002366293],[7,6.5,0.002228867],[7.5,6.5,0.002038327],[8,6.5,0.00175832],[8.5,6.5,0.001442761],[9,6.5,0.00120348],[9.5,6.5,0.000960235],[10,6.5,0.0007295162],[10.5,6.5,0.0005418395],[11,6.5,0.0004117128],[11.5,6.5,0.0003110916],[12,6.5,0.0002186289],[12.5,6.5,0.000142793],[13,6.5,9.217743e-005],[13.5,6.5,6.268011e-005],[14,6.5,4.21289e-005],[14.5,6.5,2.212571e-005],[15,6.5,9.506848e-006],[-5,7,1.014732e-005],[-4.5,7,1.697972e-005],[-4,7,2.470876e-005],[-3.5,7,3.659359e-005],[-3,7,5.908607e-005],[-2.5,7,8.028452e-005],[-2,7,0.000106393],[-1.5,7,0.0001531847],[-1,7,0.0002287426],[-0.5,7,0.0002966488],[0,7,0.0004079498],[0.5,7,0.0005449046],[1,7,0.0006918609],[1.5,7,0.0008476598],[2,7,0.001028329],[2.5,7,0.00117028],[3,7,0.001276563],[3.5,7,0.001415556],[4,7,0.001490901],[4.5,7,0.001514479],[5,7,0.00153944],[5.5,7,0.001559063],[6,7,0.001504599],[6.5,7,0.001474975],[7,7,0.001400574],[7.5,7,0.001274775],[8,7,0.00110168],[8.5,7,0.0009093013],[9,7,0.0007424355],[9.5,7,0.0006016091],[10,7,0.0004594067],[10.5,7,0.0003443522],[11,7,0.0002631525],[11.5,7,0.0001992547],[12,7,0.0001419143],[12.5,7,9.417944e-005],[13,7,6.011386e-005],[13.5,7,4.031285e-005],[14,7,2.728503e-005],[14.5,7,1.36536e-005],[15,7,5.581337e-006],[-5,7.5,4.773432e-006],[-4.5,7.5,8.189544e-006],[-4,7.5,1.229094e-005],[-3.5,7.5,1.950959e-005],[-3,7.5,2.931298e-005],[-2.5,7.5,3.960803e-005],[-2,7.5,5.622156e-005],[-1.5,7.5,7.370241e-005],[-1,7.5,0.0001152187],[-0.5,7.5,0.0001619019],[0,7.5,0.0002339767],[0.5,7.5,0.0003138408],[1,7.5,0.0003890308],[1.5,7.5,0.0004658102],[2,7.5,0.000559575],[2.5,7.5,0.0006334299],[3,7.5,0.0007044766],[3.5,7.5,0.0007791389],[4,7.5,0.0008390944],[4.5,7.5,0.00084208],[5,7.5,0.0008594495],[5.5,7.5,0.0008934814],[6,7.5,0.0008675172],[6.5,7.5,0.0008488561],[7,7.5,0.0008328978],[7.5,7.5,0.0007694918],[8,7.5,0.0006560551],[8.5,7.5,0.0005456312],[9,7.5,0.0004577678],[9.5,7.5,0.0003789428],[10,7.5,0.0002777098],[10.5,7.5,0.0002154127],[11,7.5,0.00016913],[11.5,7.5,0.0001273782],[12,7.5,9.110388e-005],[12.5,7.5,6.14494e-005],[13,7.5,3.854077e-005],[13.5,7.5,2.350409e-005],[14,7.5,1.427351e-005],[14.5,7.5,7.570889e-006],[15,7.5,3.634888e-006],[-5,8,2.058569e-006],[-4.5,8,3.739088e-006],[-4,8,6.330346e-006],[-3.5,8,2.355903e-005],[-3,8,2.449824e-005],[-2.5,8,1.696285e-005],[-2,8,2.192415e-005],[-1.5,8,3.159442e-005],[-1,8,5.137874e-005],[-0.5,8,8.260823e-005],[0,8,0.0001290211],[0.5,8,0.0001831761],[1,8,0.0002126513],[1.5,8,0.0002403995],[2,8,0.0002782417],[2.5,8,0.0003205301],[3,8,0.0003611416],[3.5,8,0.0004002437],[4,8,0.000433357],[4.5,8,0.0004303379],[5,8,0.0004315114],[5.5,8,0.0004479643],[6,8,0.0004556575],[6.5,8,0.0004624506],[7,8,0.0004656567],[7.5,8,0.0004387536],[8,8,0.0003793194],[8.5,8,0.0003158119],[9,8,0.000273984],[9.5,8,0.000230187],[10,8,0.0001620605],[10.5,8,0.0001278126],[11,8,0.0001028206],[11.5,8,7.830204e-005],[12,8,5.660197e-005],[12.5,8,3.838953e-005],[13,8,2.380503e-005],[13.5,8,1.383482e-005],[14,8,7.978304e-006],[14.5,8,4.587591e-006],[15,8,2.524582e-006],[-5,8.5,9.657881e-007],[-4.5,8.5,1.95028e-006],[-4,8.5,3.504206e-006],[-3.5,8.5,1.115178e-005],[-3,8.5,1.237526e-005],[-2.5,8.5,9.83336e-006],[-2,8.5,1.13549e-005],[-1.5,8.5,1.580112e-005],[-1,8.5,2.551574e-005],[-0.5,8.5,3.92002e-005],[0,8.5,5.750337e-005],[0.5,8.5,8.228553e-005],[1,8.5,9.91564e-005],[1.5,8.5,0.0001138667],[2,8.5,0.0001307375],[2.5,8.5,0.0001503586],[3,8.5,0.0001754117],[3.5,8.5,0.0002063177],[4,8.5,0.0002176832],[4.5,8.5,0.0002035018],[5,8.5,0.000207976],[5.5,8.5,0.0002197901],[6,8.5,0.0002308824],[6.5,8.5,0.0002466464],[7,8.5,0.0002528225],[7.5,8.5,0.000237846],[8,8.5,0.0002083612],[8.5,8.5,0.0001782277],[9,8.5,0.0001418775],[9.5,8.5,0.0001100423],[10,8.5,8.653668e-005],[10.5,8.5,7.264077e-005],[11,8.5,5.856866e-005],[11.5,8.5,4.478147e-005],[12,8.5,3.252137e-005],[12.5,8.5,2.184672e-005],[13,8.5,1.339915e-005],[13.5,8.5,7.777622e-006],[14,8.5,4.573152e-006],[14.5,8.5,2.770881e-006],[15,8.5,1.615276e-006],[-5,9,5.515082e-007],[-4.5,9,1.219677e-006],[-4,9,2.156741e-006],[-3.5,9,3.363031e-006],[-3,9,4.709053e-006],[-2.5,9,5.690698e-006],[-2,9,6.253124e-006],[-1.5,9,8.109758e-006],[-1,9,1.236813e-005],[-0.5,9,1.797675e-005],[0,9,2.47704e-005],[0.5,9,3.346485e-005],[1,9,4.255791e-005],[1.5,9,5.08377e-005],[2,9,5.890541e-005],[2.5,9,6.777852e-005],[3,9,7.833792e-005],[3.5,9,8.979053e-005],[4,9,9.45669e-005],[4.5,9,9.313145e-005],[5,9,9.711491e-005],[5.5,9,0.0001040212],[6,9,0.0001109576],[6.5,9,0.0001174021],[7,9,0.0001215219],[7.5,9,0.0001182144],[8,9,0.0001043074],[8.5,9,9.201276e-005],[9,9,7.715406e-005],[9.5,9,5.848875e-005],[10,9,4.557908e-005],[10.5,9,3.918552e-005],[11,9,3.075648e-005],[11.5,9,2.30816e-005],[12,9,1.659559e-005],[12.5,9,1.098071e-005],[13,9,6.682769e-006],[13.5,9,3.901963e-006],[14,9,2.334602e-006],[14.5,9,1.45001e-006],[15,9,8.658412e-007],[-5,9.5,3.248671e-007],[-4.5,9.5,7.438508e-007],[-4,9.5,1.299193e-006],[-3.5,9.5,1.850239e-006],[-3,9.5,2.348328e-006],[-2.5,9.5,2.745881e-006],[-2,9.5,3.154707e-006],[-1.5,9.5,4.122928e-006],[-1,9.5,5.865323e-006],[-0.5,9.5,8.040028e-006],[0,9.5,1.072878e-005],[0.5,9.5,1.430755e-005],[1,9.5,1.857546e-005],[1.5,9.5,2.299166e-005],[2,9.5,2.742243e-005],[2.5,9.5,3.187871e-005],[3,9.5,3.619584e-005],[3.5,9.5,4.015978e-005],[4,9.5,4.278464e-005],[4.5,9.5,4.390961e-005],[5,9.5,4.533535e-005],[5.5,9.5,4.775974e-005],[6,9.5,5.040043e-005],[6.5,9.5,5.233904e-005],[7,9.5,5.269405e-005],[7.5,9.5,5.100618e-005],[8,9.5,4.718537e-005],[8.5,9.5,4.190069e-005],[9,9.5,3.534603e-005],[9.5,9.5,2.71142e-005],[10,9.5,2.127897e-005],[10.5,9.5,1.77393e-005],[11,9.5,1.37143e-005],[11.5,9.5,1.022677e-005],[12,9.5,7.299217e-006],[12.5,9.5,4.809277e-006],[13,9.5,2.932456e-006],[13.5,9.5,1.712689e-006],[14,9.5,1.017901e-006],[14.5,9.5,6.285792e-007],[15,9.5,3.756647e-007],[-5,10,1.59557e-007],[-4.5,10,3.683674e-007],[-4,10,6.345202e-007],[-3.5,10,8.576232e-007],[-3,10,1.015898e-006],[-2.5,10,1.21775e-006],[-2,10,1.606562e-006],[-1.5,10,2.253995e-006],[-1,10,3.133549e-006],[-0.5,10,4.269119e-006],[0,10,5.647953e-006],[0.5,10,7.071794e-006],[1,10,8.610117e-006],[1.5,10,1.058954e-005],[2,10,1.301801e-005],[2.5,10,1.550133e-005],[3,10,1.771732e-005],[3.5,10,1.971349e-005],[4,10,2.111592e-005],[4.5,10,2.148963e-005],[5,10,2.143981e-005],[5.5,10,2.149826e-005],[6,10,2.15542e-005],[6.5,10,2.161197e-005],[7,10,2.154276e-005],[7.5,10,2.125465e-005],[8,10,2.065819e-005],[8.5,10,1.906993e-005],[9,10,1.598837e-005],[9.5,10,1.210061e-005],[10,10,8.911344e-006],[10.5,10,6.755228e-006],[11,10,5.117238e-006],[11.5,10,3.8145e-006],[12,10,2.716026e-006],[12.5,10,1.802211e-006],[13,10,1.108003e-006],[13.5,10,6.421545e-007],[14,10,3.708128e-007],[14.5,10,2.218171e-007],[15,10,1.302674e-007],[-5,10.5,5.939211e-008],[-4.5,10.5,1.375441e-007],[-4,10.5,2.368123e-007],[-3.5,10.5,3.199517e-007],[-3,10.5,3.918598e-007],[-2.5,10.5,5.319167e-007],[-2,10.5,8.162968e-007],[-1.5,10.5,1.252473e-006],[-1,10.5,1.920157e-006],[-0.5,10.5,3.019602e-006],[0,10.5,4.260673e-006],[0.5,10.5,4.789419e-006],[1,10.5,4.645738e-006],[1.5,10.5,4.876109e-006],[2,10.5,5.893107e-006],[2.5,10.5,7.186595e-006],[3,10.5,8.285914e-006],[3.5,10.5,9.214796e-006],[4,10.5,9.905516e-006],[4.5,10.5,1.006634e-005],[5,10.5,9.797952e-006],[5.5,10.5,9.297572e-006],[6,10.5,8.604366e-006],[6.5,10.5,8.006683e-006],[7,10.5,7.795449e-006],[7.5,10.5,7.979307e-006],[8,10.5,8.314752e-006],[8.5,10.5,8.226698e-006],[9,10.5,7.191361e-006],[9.5,10.5,5.407573e-006],[10,10.5,3.653434e-006],[10.5,10.5,2.429904e-006],[11,10.5,1.684655e-006],[11.5,10.5,1.198666e-006],[12,10.5,8.391078e-007],[12.5,10.5,5.598685e-007],[13,10.5,3.468033e-007],[13.5,10.5,1.986883e-007],[14,10.5,1.10408e-007],[14.5,10.5,6.301073e-008],[15,10.5,3.589721e-008],[-5,11,1.626954e-008],[-4.5,11,3.784139e-008],[-4,11,6.6018e-008],[-3.5,11,9.33614e-008],[-3,11,1.290908e-007],[-2.5,11,2.078178e-007],[-2,11,3.601706e-007],[-1.5,11,6.170417e-007],[-1,11,1.145169e-006],[-0.5,11,2.184055e-006],[0,11,3.335304e-006],[0.5,11,3.554577e-006],[1,11,2.809516e-006],[1.5,11,2.262392e-006],[2,11,2.466164e-006],[2.5,11,3.00369e-006],[3,11,3.435902e-006],[3.5,11,3.733008e-006],[4,11,3.958295e-006],[4.5,11,4.030607e-006],[5,11,3.907656e-006],[5.5,11,3.580934e-006],[6,11,3.074849e-006],[6.5,11,2.619389e-006],[7,11,2.467472e-006],[7.5,11,2.65895e-006],[8,11,3.030965e-006],[8.5,11,3.249551e-006],[9,11,2.989995e-006],[9.5,11,2.270085e-006],[10,11,1.450295e-006],[10.5,11,8.450239e-007],[11,11,5.000279e-007],[11.5,11,3.171134e-007],[12,11,2.10882e-007],[12.5,11,1.394817e-007],[13,11,8.662414e-008],[13.5,11,4.909607e-008],[14,11,2.629811e-008],[14.5,11,1.426516e-008],[15,11,7.831502e-009],[-5,11.5,3.254276e-009],[-4.5,11.5,7.630869e-009],[-4,11.5,1.369958e-008],[-3.5,11.5,2.116599e-008],[-3,11.5,3.486278e-008],[-2.5,11.5,6.614921e-008],[-2,11.5,1.261958e-007],[-1.5,11.5,2.443782e-007],[-1,11.5,5.46815e-007],[-0.5,11.5,1.189942e-006],[0,11.5,1.903684e-006],[0.5,11.5,2.002463e-006],[1,11.5,1.44495e-006],[1.5,11.5,9.607321e-007],[2,11.5,9.170597e-007],[2.5,11.5,1.084515e-006],[3,11.5,1.204591e-006],[3.5,11.5,1.243813e-006],[4,11.5,1.26525e-006],[4.5,11.5,1.277832e-006],[5,11.5,1.245216e-006],[5.5,11.5,1.122563e-006],[6,11.5,9.139102e-007],[6.5,11.5,7.261132e-007],[7,11.5,6.764284e-007],[7.5,11.5,7.920044e-007],[8,11.5,1.006558e-006],[8.5,11.5,1.172341e-006],[9,11.5,1.13404e-006],[9.5,11.5,8.754691e-007],[10,11.5,5.413918e-007],[10.5,11.5,2.820651e-007],[11,11.5,1.377824e-007],[11.5,11.5,7.200651e-008],[12,11.5,4.275597e-008],[12.5,11.5,2.725857e-008],[13,11.5,1.683747e-008],[13.5,11.5,9.461957e-009],[14,11.5,4.921004e-009],[14.5,11.5,2.549141e-009],[15,11.5,1.348085e-009],[-5,12,4.74764e-010],[-4.5,12,1.130408e-009],[-4,12,2.137466e-009],[-3.5,12,3.785364e-009],[-3,12,7.598643e-009],[-2.5,12,1.645695e-008],[-2,12,3.368567e-008],[-1.5,12,7.209029e-008],[-1,12,1.824898e-007],[-0.5,12,4.248474e-007],[0,12,6.948243e-007],[0.5,12,7.304255e-007],[1,12,5.13449e-007],[1.5,12,3.179007e-007],[2,12,2.831286e-007],[2.5,12,3.25199e-007],[3,12,3.477023e-007],[3.5,12,3.363559e-007],[4,12,3.202634e-007],[4.5,12,3.140188e-007],[5,12,3.047227e-007],[5.5,12,2.710103e-007],[6,12,2.128583e-007],[6.5,12,1.633735e-007],[7,12,1.589919e-007],[7.5,12,2.113697e-007],[8,12,3.014538e-007],[8.5,12,3.773804e-007],[9,12,3.799712e-007],[9.5,12,2.980615e-007],[10,12,1.815911e-007],[10.5,12,8.801783e-008],[11,12,3.640058e-008],[11.5,12,1.485838e-008],[12,12,7.154605e-009],[12.5,12,4.158986e-009],[13,12,2.51441e-009],[13.5,12,1.401678e-009],[14,12,7.136911e-010],[14.5,12,3.564351e-010],[15,12,1.824956e-010],[-5,12.5,5.061484e-011],[-4.5,12.5,1.240092e-010],[-4,12.5,2.561736e-010],[-3.5,12.5,5.454587e-010],[-3,12.5,1.32032e-009],[-2.5,12.5,3.133892e-009],[-2,12.5,6.673219e-009],[-1.5,12.5,1.496287e-008],[-1,12.5,3.977527e-008],[-0.5,12.5,9.483966e-008],[0,12.5,1.564616e-007],[0.5,12.5,1.652857e-007],[1,12.5,1.174512e-007],[1.5,12.5,7.51618e-008],[2,12.5,6.924778e-008],[2.5,12.5,7.915975e-008],[3,12.5,8.183721e-008],[3.5,12.5,7.434197e-008],[4,12.5,6.542432e-008],[4.5,12.5,6.060696e-008],[5,12.5,5.713078e-008],[5.5,12.5,4.956183e-008],[6,12.5,3.78953e-008],[6.5,12.5,2.932162e-008],[7,12.5,3.204435e-008],[7.5,12.5,4.97383e-008],[8,12.5,7.830426e-008],[8.5,12.5,1.032128e-007],[9,12.5,1.06714e-007],[9.5,12.5,8.467899e-008],[10,12.5,5.134202e-008],[10.5,12.5,2.400681e-008],[11,12.5,8.956881e-009],[11.5,12.5,2.939062e-009],[12,12.5,1.054272e-009],[12.5,12.5,5.05023e-010],[13,12.5,2.878134e-010],[13.5,12.5,1.583606e-010],[14,12.5,7.949094e-011],[14.5,12.5,3.872393e-011],[15,12.5,1.937346e-011],[-5,13,3.963356e-012],[-4.5,13,1.023587e-011],[-4,13,2.430987e-011],[-3.5,13,6.400965e-011],[-3,13,1.798711e-010],[-2.5,13,4.509524e-010],[-2,13,9.676819e-010],[-1.5,13,2.13007e-009],[-1,13,5.534663e-009],[-0.5,13,1.307451e-008],[0,13,2.15861e-008],[0.5,13,2.311096e-008],[1,13,1.732426e-008],[1.5,13,1.279205e-008],[2,13,1.338282e-008],[2.5,13,1.557581e-008],[3,13,1.569281e-008],[3.5,13,1.352254e-008],[4,13,1.098771e-008],[4.5,13,9.38247e-009],[5,13,8.330693e-009],[5.5,13,6.924195e-009],[6,13,5.17951e-009],[6.5,13,4.240384e-009],[7,13,5.532268e-009],[7.5,13,9.937467e-009],[8,13,1.676941e-008],[8.5,13,2.281558e-008],[9,13,2.395821e-008],[9.5,13,1.914639e-008],[10,13,1.160017e-008],[10.5,13,5.343418e-009],[11,13,1.897689e-009],[11.5,13,5.447933e-010],[12,13,1.474408e-010],[12.5,13,5.153118e-011],[13,13,2.551202e-011],[13.5,13,1.360897e-011],[14,13,6.762243e-012],[14.5,13,3.253054e-012],[15,13,1.609715e-012],[-5,13.5,2.30342e-013],[-4.5,13.5,6.529597e-013],[-4,13.5,1.88248e-012],[-3.5,13.5,6.072555e-012],[-3,13.5,1.889188e-011],[-2.5,13.5,4.871282e-011],[-2,13.5,1.028721e-010],[-1.5,13.5,2.09902e-010],[-1,13.5,4.957986e-010],[-0.5,13.5,1.119515e-009],[0,13.5,1.84072e-009],[0.5,13.5,2.044551e-009],[1,13.5,1.755046e-009],[1.5,13.5,1.692653e-009],[2,13.5,2.091021e-009],[2.5,13.5,2.487593e-009],[3,13.5,2.455133e-009],[3.5,13.5,2.025635e-009],[4,13.5,1.526169e-009],[4.5,13.5,1.182029e-009],[5,13.5,9.629759e-010],[5.5,13.5,7.544537e-010],[6,13.5,5.556205e-010],[6.5,13.5,5.052242e-010],[7,13.5,8.04138e-010],[7.5,13.5,1.61398e-009],[8,13.5,2.842583e-009],[8.5,13.5,3.937922e-009],[9,13.5,4.170867e-009],[9.5,13.5,3.346788e-009],[10,13.5,2.028567e-009],[10.5,13.5,9.291936e-010],[11,13.5,3.233228e-010],[11.5,13.5,8.71463e-011],[12,13.5,1.964817e-011],[12.5,13.5,4.809219e-012],[13,13.5,1.808829e-012],[13.5,13.5,8.925797e-013],[14,13.5,4.383798e-013],[14.5,13.5,2.108244e-013],[15,13.5,1.046368e-013],[-5,14,1.013541e-014],[-4.5,14,3.340682e-014],[-4,14,1.2064e-013],[-3.5,14,4.565489e-013],[-3,14,1.509866e-012],[-2.5,14,3.942094e-012],[-2,14,8.116482e-012],[-1.5,14,1.489092e-011],[-1,14,2.983633e-011],[-0.5,14,6.14068e-011],[0,14,1.009455e-010],[0.5,14,1.246442e-010],[1,14,1.420251e-010],[1.5,14,1.909602e-010],[2,14,2.685317e-010],[2.5,14,3.231137e-010],[3,14,3.133216e-010],[3.5,14,2.492938e-010],[4,14,1.751883e-010],[4.5,14,1.22151e-010],[5,14,8.967145e-011],[5.5,14,6.54986e-011],[6,14,4.796917e-011],[6.5,14,5.019952e-011],[7,14,9.542501e-011],[7.5,14,2.060282e-010],[8,14,3.719861e-010],[8.5,14,5.2056e-010],[9,14,5.539841e-010],[9.5,14,4.455552e-010],[10,14,2.702147e-010],[10.5,14,1.235321e-010],[11,14,4.264583e-011],[11.5,14,1.120039e-011],[12,14,2.310756e-012],[12.5,14,4.323177e-013],[13,14,1.097764e-013],[13.5,14,4.53347e-014],[14,14,2.168748e-014],[14.5,14,1.05465e-014],[15,14,5.328708e-015],[-5,14.5,3.49186e-016],[-4.5,14.5,1.419874e-015],[-4,14.5,6.330915e-015],[-3.5,14.5,2.661892e-014],[-3,14.5,9.102899e-014],[-2.5,14.5,2.38888e-013],[-2,14.5,4.812969e-013],[-1.5,14.5,8.009783e-013],[-1,14.5,1.317984e-012],[-0.5,14.5,2.36639e-012],[0,14.5,4.069408e-012],[0.5,14.5,6.517696e-012],[1,14.5,1.100025e-011],[1.5,14.5,1.887314e-011],[2,14.5,2.82244e-011],[2.5,14.5,3.398557e-011],[3,14.5,3.248209e-011],[3.5,14.5,2.509547e-011],[4,14.5,1.658606e-011],[4.5,14.5,1.041237e-011],[5,14.5,6.801406e-012],[5.5,14.5,4.590458e-012],[6,14.5,3.37787e-012],[6.5,14.5,4.118361e-012],[7,14.5,8.967022e-012],[7.5,14.5,2.024845e-011],[8,14.5,3.708719e-011],[8.5,14.5,5.22005e-011],[9,14.5,5.570255e-011],[9.5,14.5,4.485968e-011],[10,14.5,2.721835e-011],[10.5,14.5,1.24356e-011],[11,14.5,4.280353e-012],[11.5,14.5,1.11297e-012],[12,14.5,2.213193e-013],[12.5,14.5,3.58334e-014],[13,14.5,6.220186e-015],[13.5,14.5,1.851099e-015],[14,14.5,8.245877e-016],[14.5,14.5,4.087833e-016],[15,14.5,2.133126e-016],[-5,15,9.850298e-018],[-4.5,15,5.074215e-017],[-4,15,2.650157e-016],[-3.5,15,1.18369e-015],[-3,15,4.118191e-015],[-2.5,15,1.083326e-014],[-2,15,2.156172e-014],[-1.5,15,3.382703e-014],[-1,15,4.819355e-014],[-0.5,15,7.921199e-014],[0,15,1.632077e-013],[0.5,15,3.748396e-013],[1,15,8.297885e-013],[1.5,15,1.571291e-012],[2,15,2.391722e-012],[2.5,15,2.869854e-012],[3,15,2.71391e-012],[3.5,15,2.050923e-012],[4,15,1.289397e-012],[4.5,15,7.347768e-013],[5,15,4.235302e-013],[5.5,15,2.607333e-013],[6,15,1.934308e-013],[6.5,15,2.723947e-013],[7,15,6.528505e-013],[7.5,15,1.515129e-012],[8,15,2.798858e-012],[8.5,15,3.952894e-012],[9,15,4.224813e-012],[9.5,15,3.405107e-012],[10,15,2.066694e-012],[10.5,15,9.440943e-013],[11,15,3.246085e-013],[11.5,15,8.40855e-014],[12,15,1.648469e-014],[12.5,15,2.506434e-015],[13,15,3.385395e-016],[13.5,15,6.531208e-017],[14,15,2.451912e-017],[14.5,15,1.236827e-017],[15,15,6.746042e-018]],"colors":[[0,0,1,1]],"centers":[[-5,-5,1.022358e-005],[-4.5,-5,1.103743e-005],[-4,-5,1.318789e-005],[-3.5,-5,1.733092e-005],[-3,-5,2.381306e-005],[-2.5,-5,3.314449e-005],[-2,-5,4.359418e-005],[-1.5,-5,5.19521e-005],[-1,-5,6.044171e-005],[-0.5,-5,7.313126e-005],[0,-5,8.640253e-005],[0.5,-5,0.0001005138],[1,-5,0.0001102444],[1.5,-5,0.0001154121],[2,-5,0.0001190983],[2.5,-5,0.0001233736],[3,-5,0.0001256147],[3.5,-5,0.00012259],[4,-5,0.0001157894],[4.5,-5,0.0001086386],[5,-5,0.0001049353],[5.5,-5,8.828189e-005],[6,-5,7.543493e-005],[6.5,-5,6.864825e-005],[7,-5,5.840666e-005],[7.5,-5,4.45579e-005],[8,-5,3.183114e-005],[8.5,-5,2.234906e-005],[9,-5,1.559137e-005],[9.5,-5,1.020451e-005],[10,-5,5.836081e-006],[10.5,-5,2.841148e-006],[11,-5,1.225685e-006],[11.5,-5,6.253013e-007],[12,-5,5.327451e-007],[12.5,-5,3.598875e-007],[13,-5,1.156378e-007],[13.5,-5,2.102437e-008],[14,-5,5.226523e-009],[14.5,-5,3.040046e-009],[15,-5,2.85077e-009],[-5,-4.5,2.069899e-005],[-4.5,-4.5,2.323213e-005],[-4,-4.5,2.796962e-005],[-3.5,-4.5,3.652863e-005],[-3,-4.5,4.949964e-005],[-2.5,-4.5,6.759587e-005],[-2,-4.5,8.810722e-005],[-1.5,-4.5,0.0001060953],[-1,-4.5,0.0001253808],[-0.5,-4.5,0.0001508964],[0,-4.5,0.0001732802],[0.5,-4.5,0.0001988866],[1,-4.5,0.0002193458],[1.5,-4.5,0.0002361014],[2,-4.5,0.0002513785],[2.5,-4.5,0.000254854],[3,-4.5,0.0002701528],[3.5,-4.5,0.0002624239],[4,-4.5,0.0002550682],[4.5,-4.5,0.000227719],[5,-4.5,0.0002172898],[5.5,-4.5,0.0001897597],[6,-4.5,0.0001542671],[6.5,-4.5,0.0001306565],[7,-4.5,0.0001076691],[7.5,-4.5,8.379408e-005],[8,-4.5,6.203192e-005],[8.5,-4.5,4.501045e-005],[9,-4.5,3.251125e-005],[9.5,-4.5,2.241555e-005],[10,-4.5,1.39816e-005],[10.5,-4.5,7.73611e-006],[11,-4.5,3.969961e-006],[11.5,-4.5,2.757731e-006],[12,-4.5,3.067434e-006],[12.5,-4.5,2.218886e-006],[13,-4.5,7.071512e-007],[13.5,-4.5,1.201952e-007],[14,-4.5,2.868395e-008],[14.5,-4.5,1.734266e-008],[15,-4.5,1.631938e-008],[-5,-4,3.352584e-005],[-4.5,-4,4.049542e-005],[-4,-4,5.214016e-005],[-3.5,-4,7.164631e-005],[-3,-4,9.75771e-005],[-2.5,-4,0.0001290569],[-2,-4,0.0001648147],[-1.5,-4,0.0002035393],[-1,-4,0.0002423821],[-0.5,-4,0.000277658],[0,-4,0.0003232002],[0.5,-4,0.0003687174],[1,-4,0.0004159762],[1.5,-4,0.0004651498],[2,-4,0.0004928176],[2.5,-4,0.0004928212],[3,-4,0.0005251947],[3.5,-4,0.0005273412],[4,-4,0.0005226753],[4.5,-4,0.0004599563],[5,-4,0.0004064105],[5.5,-4,0.0003559696],[6,-4,0.0003009211],[6.5,-4,0.0002474041],[7,-4,0.0002011052],[7.5,-4,0.0001592811],[8,-4,0.000121623],[8.5,-4,9.098976e-005],[9,-4,6.744826e-005],[9.5,-4,4.823787e-005],[10,-4,3.217364e-005],[10.5,-4,1.960548e-005],[11,-4,1.114271e-005],[11.5,-4,7.964444e-006],[12,-4,8.288611e-006],[12.5,-4,5.854569e-006],[13,-4,1.943029e-006],[13.5,-4,3.967965e-007],[14,-4,1.23337e-007],[14.5,-4,7.886126e-008],[15,-4,7.274446e-008],[-5,-3.5,4.694144e-005],[-4.5,-3.5,6.258528e-005],[-4,-3.5,9.101545e-005],[-3.5,-3.5,0.00013642],[-3,-3.5,0.0001868602],[-2.5,-3.5,0.0002384381],[-2,-3.5,0.0002940826],[-1.5,-3.5,0.0003653084],[-1,-3.5,0.0004399769],[-0.5,-3.5,0.0005080014],[0,-3.5,0.0005885305],[0.5,-3.5,0.0006682887],[1,-3.5,0.0007498874],[1.5,-3.5,0.000833279],[2,-3.5,0.0008924222],[2.5,-3.5,0.0009123718],[3,-3.5,0.0009428606],[3.5,-3.5,0.0009573893],[4,-3.5,0.0009415169],[4.5,-3.5,0.0008548081],[5,-3.5,0.0007581883],[5.5,-3.5,0.0006717056],[6,-3.5,0.0005796183],[6.5,-3.5,0.0004669309],[7,-3.5,0.0003775536],[7.5,-3.5,0.0003008256],[8,-3.5,0.0002325377],[8.5,-3.5,0.0001762972],[9,-3.5,0.0001313529],[9.5,-3.5,9.535612e-005],[10,-3.5,6.698941e-005],[10.5,-3.5,4.461935e-005],[11,-3.5,2.761806e-005],[11.5,-3.5,1.714831e-005],[12,-3.5,1.286095e-005],[12.5,-3.5,8.092175e-006],[13,-3.5,3.111339e-006],[13.5,-3.5,9.944935e-007],[14,-3.5,4.363978e-007],[14.5,-3.5,2.895069e-007],[15,-3.5,2.559674e-007],[-5,-3,6.283291e-005],[-4.5,-3,9.251746e-005],[-4,-3,0.0001471142],[-3.5,-3,0.0002297969],[-3,-3,0.0003121307],[-2.5,-3,0.0004001023],[-2,-3,0.0005044413],[-1.5,-3,0.0006330485],[-1,-3,0.0007576831],[-0.5,-3,0.0008751507],[0,-3,0.001006732],[0.5,-3,0.001166384],[1,-3,0.001326532],[1.5,-3,0.001486583],[2,-3,0.001629656],[2.5,-3,0.001650302],[3,-3,0.001670886],[3.5,-3,0.001664893],[4,-3,0.001608517],[4.5,-3,0.001487143],[5,-3,0.001343528],[5.5,-3,0.001183867],[6,-3,0.001019062],[6.5,-3,0.0008525469],[7,-3,0.0006903602],[7.5,-3,0.0005486452],[8,-3,0.0004263885],[8.5,-3,0.0003392127],[9,-3,0.0002492083],[9.5,-3,0.0001753216],[10,-3,0.0001275737],[10.5,-3,9.32021e-005],[11,-3,6.687154e-005],[11.5,-3,3.787612e-005],[12,-3,2.088358e-005],[12.5,-3,1.170257e-005],[13,-3,5.556998e-006],[13.5,-3,2.539993e-006],[14,-3,1.321974e-006],[14.5,-3,8.805005e-007],[15,-3,7.269303e-007],[-5,-2.5,8.317419e-005],[-4.5,-2.5,0.0001287157],[-4,-2.5,0.0002019633],[-3.5,-2.5,0.0003048374],[-3,-2.5,0.0004271723],[-2.5,-2.5,0.0005840205],[-2,-2.5,0.0007828551],[-1.5,-2.5,0.001003632],[-1,-2.5,0.0012005],[-0.5,-2.5,0.001396805],[0,-2.5,0.001652539],[0.5,-2.5,0.001935821],[1,-2.5,0.002161358],[1.5,-2.5,0.002413663],[2,-2.5,0.002598478],[2.5,-2.5,0.002700389],[3,-2.5,0.002729656],[3.5,-2.5,0.002745303],[4,-2.5,0.002672856],[4.5,-2.5,0.002461533],[5,-2.5,0.002237034],[5.5,-2.5,0.001982393],[6,-2.5,0.001717946],[6.5,-2.5,0.001445418],[7,-2.5,0.001193795],[7.5,-2.5,0.0009578428],[8,-2.5,0.0007434991],[8.5,-2.5,0.0005971643],[9,-2.5,0.0004381583],[9.5,-2.5,0.0003055494],[10,-2.5,0.0002236774],[10.5,-2.5,0.0001677915],[11,-2.5,0.0001164451],[11.5,-2.5,6.919775e-005],[12,-2.5,4.006087e-005],[12.5,-2.5,2.291202e-005],[13,-2.5,1.204298e-005],[13.5,-2.5,6.205087e-006],[14,-2.5,3.464526e-006],[14.5,-2.5,2.271179e-006],[15,-2.5,1.70555e-006],[-5,-2,0.0001093828],[-4.5,-2,0.0001718271],[-4,-2,0.0002713594],[-3.5,-2,0.0004102671],[-3,-2,0.000606637],[-2.5,-2,0.0008770715],[-2,-2,0.001139808],[-1.5,-2,0.001483015],[-1,-2,0.001788957],[-0.5,-2,0.002087827],[0,-2,0.002536748],[0.5,-2,0.0029928],[1,-2,0.003341079],[1.5,-2,0.003675518],[2,-2,0.003913896],[2.5,-2,0.004107998],[3,-2,0.004225537],[3.5,-2,0.004248867],[4,-2,0.004191169],[4.5,-2,0.003934869],[5,-2,0.003611415],[5.5,-2,0.003176386],[6,-2,0.002791955],[6.5,-2,0.002387549],[7,-2,0.001954116],[7.5,-2,0.001578444],[8,-2,0.001234379],[8.5,-2,0.0009568912],[9,-2,0.0007223485],[9.5,-2,0.0005203434],[10,-2,0.0003706378],[10.5,-2,0.0002685387],[11,-2,0.0001847422],[11.5,-2,0.0001216623],[12,-2,7.370643e-005],[12.5,-2,4.366912e-005],[13,-2,2.452365e-005],[13.5,-2,1.356168e-005],[14,-2,7.896796e-006],[14.5,-2,4.996446e-006],[15,-2,3.35533e-006],[-5,-1.5,0.0001523374],[-4.5,-1.5,0.0002392302],[-4,-1.5,0.0003742393],[-3.5,-1.5,0.0005713367],[-3,-1.5,0.0008743363],[-2.5,-1.5,0.001272012],[-2,-1.5,0.001629646],[-1.5,-1.5,0.002085713],[-1,-1.5,0.002547823],[-0.5,-1.5,0.003040961],[0,-1.5,0.003670539],[0.5,-1.5,0.004358381],[1,-1.5,0.004943821],[1.5,-1.5,0.00544283],[2,-1.5,0.00580303],[2.5,-1.5,0.006083431],[3,-1.5,0.006255389],[3.5,-1.5,0.006263445],[4,-1.5,0.006231996],[4.5,-1.5,0.005965018],[5,-1.5,0.005470508],[5.5,-1.5,0.004870108],[6,-1.5,0.004405541],[6.5,-1.5,0.003792814],[7,-1.5,0.003072205],[7.5,-1.5,0.002482503],[8,-1.5,0.001972175],[8.5,-1.5,0.001516801],[9,-1.5,0.001173589],[9.5,-1.5,0.0008787551],[10,-1.5,0.0006026304],[10.5,-1.5,0.0004110658],[11,-1.5,0.0002887338],[11.5,-1.5,0.0002069897],[12,-1.5,0.0001298236],[12.5,-1.5,7.481874e-005],[13,-1.5,4.448937e-005],[13.5,-1.5,2.616374e-005],[14,-1.5,1.557472e-005],[14.5,-1.5,9.368237e-006],[15,-1.5,5.62355e-006],[-5,-1,0.0002136888],[-4.5,-1,0.0003333431],[-4,-1,0.0005013258],[-3.5,-1,0.0007607086],[-3,-1,0.001157601],[-2.5,-1,0.001646307],[-2,-1,0.002114639],[-1.5,-1,0.002712109],[-1,-1,0.003365271],[-0.5,-1,0.004141956],[0,-1,0.00507151],[0.5,-1,0.006020527],[1,-1,0.006887507],[1.5,-1,0.007648244],[2,-1,0.008223148],[2.5,-1,0.008626256],[3,-1,0.008855464],[3.5,-1,0.008871207],[4,-1,0.00884754],[4.5,-1,0.008500398],[5,-1,0.007864795],[5.5,-1,0.007073918],[6,-1,0.006330473],[6.5,-1,0.005483005],[7,-1,0.004553836],[7.5,-1,0.003709831],[8,-1,0.002964412],[8.5,-1,0.002274409],[9,-1,0.001699988],[9.5,-1,0.001265323],[10,-1,0.0008921861],[10.5,-1,0.0006048291],[11,-1,0.0004253611],[11.5,-1,0.000295162],[12,-1,0.0001942911],[12.5,-1,0.0001176085],[13,-1,7.219736e-005],[13.5,-1,4.442909e-005],[14,-1,2.67895e-005],[14.5,-1,1.538571e-005],[15,-1,8.441748e-006],[-5,-0.5,0.0002770707],[-4.5,-0.5,0.0004275369],[-4,-0.5,0.000636605],[-3.5,-0.5,0.0009603131],[-3,-0.5,0.001425674],[-2.5,-0.5,0.001933799],[-2,-0.5,0.002548859],[-1.5,-0.5,0.003309109],[-1,-0.5,0.004211568],[-0.5,-0.5,0.005276743],[0,-0.5,0.006549423],[0.5,-0.5,0.007879692],[1,-0.5,0.009108986],[1.5,-0.5,0.01014582],[2,-0.5,0.01098782],[2.5,-0.5,0.01164978],[3,-0.5,0.01205101],[3.5,-0.5,0.01212336],[4,-0.5,0.01200153],[4.5,-0.5,0.01161393],[5,-0.5,0.0109102],[5.5,-0.5,0.009806607],[6,-0.5,0.008608765],[6.5,-0.5,0.007440581],[7,-0.5,0.006348599],[7.5,-0.5,0.00516997],[8,-0.5,0.004149278],[8.5,-0.5,0.003169118],[9,-0.5,0.002349548],[9.5,-0.5,0.001728037],[10,-0.5,0.001260719],[10.5,-0.5,0.0008905967],[11,-0.5,0.0006411481],[11.5,-0.5,0.0004184243],[12,-0.5,0.000272737],[12.5,-0.5,0.0001730736],[13,-0.5,0.0001066119],[13.5,-0.5,6.724011e-005],[14,-0.5,4.109306e-005],[14.5,-0.5,2.316816e-005],[15,-0.5,1.219354e-005],[-5,0,0.0003209224],[-4.5,0,0.0005134235],[-4,0,0.0007811526],[-3.5,0,0.001155516],[-3,0,0.001656024],[-2.5,0,0.002254551],[-2,0,0.003031094],[-1.5,0,0.003888661],[-1,0,0.005040784],[-0.5,0,0.006334898],[0,0,0.007765371],[0.5,0,0.009466299],[1,0,0.01119199],[1.5,0,0.01256562],[2,0,0.01372529],[2.5,0,0.01473876],[3,0,0.01540706],[3.5,0,0.01564884],[4,0,0.01544846],[4.5,0,0.01491693],[5,0,0.01399857],[5.5,0,0.01274496],[6,0,0.0112606],[6.5,0,0.009763407],[7,0,0.008276966],[7.5,0,0.006787405],[8,0,0.005504527],[8.5,0,0.004269429],[9,0,0.003221611],[9.5,0,0.002376397],[10,0,0.001716942],[10.5,0,0.00123333],[11,0,0.0008741664],[11.5,0,0.0005769234],[12,0,0.0003682042],[12.5,0,0.0002346456],[13,0,0.0001459708],[13.5,0,9.256292e-005],[14,0,5.688789e-005],[14.5,0,3.211615e-005],[15,0,1.696517e-005],[-5,0.5,0.0003593504],[-4.5,0.5,0.0005729354],[-4,0.5,0.0008808575],[-3.5,0.5,0.001294177],[-3,0.5,0.001845379],[-2.5,0.5,0.002540605],[-2,0.5,0.003430831],[-1.5,0.5,0.004458918],[-1,0.5,0.00577658],[-0.5,0.5,0.007278748],[0,0.5,0.00888849],[0.5,0.5,0.0108683],[1,0.5,0.01292731],[1.5,0.5,0.01464605],[2,0.5,0.01622435],[2.5,0.5,0.0176887],[3,0.5,0.01859451],[3.5,0.5,0.01897665],[4,0.5,0.01891237],[4.5,0.5,0.01822779],[5,0.5,0.01715785],[5.5,0.5,0.01579985],[6,0.5,0.01400261],[6.5,0.5,0.01221664],[7,0.5,0.01039986],[7.5,0.5,0.008624972],[8,0.5,0.006918869],[8.5,0.5,0.00535449],[9,0.5,0.004114065],[9.5,0.5,0.003094596],[10,0.5,0.002246772],[10.5,0.5,0.001618595],[11,0.5,0.001131231],[11.5,0.5,0.0007361925],[12,0.5,0.000475811],[12.5,0.5,0.0003090613],[13,0.5,0.0001923105],[13.5,0.5,0.0001255382],[14,0.5,7.590792e-005],[14.5,0.5,4.042645e-005],[15,0.5,2.15269e-005],[-5,1,0.0003643982],[-4.5,1,0.0005795374],[-4,1,0.0008874646],[-3.5,1,0.001318627],[-3,1,0.001903539],[-2.5,1,0.002640061],[-2,1,0.003648189],[-1.5,1,0.004888502],[-1,1,0.006412107],[-0.5,1,0.00805085],[0,1,0.009901511],[0.5,1,0.01209675],[1,1,0.01435261],[1.5,1,0.01647073],[2,1,0.01847414],[2.5,1,0.0200619],[3,1,0.0211264],[3.5,1,0.0216377],[4,1,0.0216244],[4.5,1,0.02106634],[5,1,0.01994348],[5.5,1,0.0184497],[6,1,0.0165868],[6.5,1,0.01455675],[7,1,0.01253109],[7.5,1,0.01031779],[8,1,0.008187769],[8.5,1,0.006428245],[9,1,0.005023892],[9.5,1,0.003807913],[10,1,0.002770599],[10.5,1,0.002007507],[11,1,0.00140066],[11.5,1,0.0009044391],[12,1,0.0005990935],[12.5,1,0.0003827609],[13,1,0.0002416837],[13.5,1,0.0001558393],[14,1,9.121562e-005],[14.5,1,4.713896e-005],[15,1,2.500576e-005],[-5,1.5,0.0003455108],[-4.5,1.5,0.0005551279],[-4,1.5,0.0008628383],[-3.5,1.5,0.001296526],[-3,1.5,0.00186034],[-2.5,1.5,0.002663968],[-2,1.5,0.00371794],[-1.5,1.5,0.005032463],[-1,1.5,0.006645836],[-0.5,1.5,0.008366022],[0,1.5,0.01037911],[0.5,1.5,0.01270083],[1,1.5,0.01515518],[1.5,1.5,0.01752976],[2,1.5,0.01964824],[2.5,1.5,0.02122781],[3,1.5,0.02255295],[3.5,1.5,0.02335948],[4,1.5,0.02360849],[4.5,1.5,0.02303916],[5,1.5,0.02184367],[5.5,1.5,0.02030874],[6,1.5,0.01852154],[6.5,1.5,0.01637102],[7,1.5,0.01390242],[7.5,1.5,0.01151755],[8,1.5,0.009351893],[8.5,1.5,0.007391411],[9,1.5,0.005814216],[9.5,1.5,0.004422875],[10,1.5,0.003235861],[10.5,1.5,0.002316762],[11,1.5,0.001612895],[11.5,1.5,0.001057941],[12,1.5,0.0007024868],[12.5,1.5,0.0004524757],[13,1.5,0.0002974331],[13.5,1.5,0.0001774191],[14,1.5,9.890705e-005],[14.5,1.5,5.325396e-005],[15,1.5,2.803842e-005],[-5,2,0.0003409614],[-4.5,2,0.0005307131],[-4,2,0.0008160742],[-3.5,2,0.001263079],[-3,2,0.0018135],[-2.5,2,0.002597569],[-2,2,0.003628985],[-1.5,2,0.004940581],[-1,2,0.006500918],[-0.5,2,0.008184075],[0,2,0.01021786],[0.5,2,0.01263407],[1,2,0.01503569],[1.5,2,0.01736504],[2,2,0.01951774],[2.5,2,0.0212557],[3,2,0.02268321],[3.5,2,0.02363331],[4,2,0.02414834],[4.5,2,0.02376101],[5,2,0.02266241],[5.5,2,0.0212651],[6,2,0.01941289],[6.5,2,0.01724573],[7,2,0.01475658],[7.5,2,0.01238161],[8,2,0.01026321],[8.5,2,0.008163632],[9,2,0.006462851],[9.5,2,0.004850509],[10,2,0.003579124],[10.5,2,0.002564469],[11,2,0.001775951],[11.5,2,0.001188788],[12,2,0.0007933318],[12.5,2,0.0005130534],[13,2,0.0003411908],[13.5,2,0.0002005631],[14,2,0.0001106596],[14.5,2,5.930436e-005],[15,2,3.103482e-005],[-5,2.5,0.0003000559],[-4.5,2.5,0.0004768276],[-4,2.5,0.0007216992],[-3.5,2.5,0.001124901],[-3,2.5,0.001645392],[-2.5,2.5,0.002346262],[-2,2.5,0.003320437],[-1.5,2.5,0.004548607],[-1,2.5,0.006013018],[-0.5,2.5,0.007680389],[0,2.5,0.009529932],[0.5,2.5,0.01184772],[1,2.5,0.01408372],[1.5,2.5,0.01629614],[2,2.5,0.01836063],[2.5,2.5,0.02033126],[3,2.5,0.02178055],[3.5,2.5,0.02278999],[4,2.5,0.02336491],[4.5,2.5,0.02317152],[5,2.5,0.02239033],[5.5,2.5,0.02120073],[6,2.5,0.01938601],[6.5,2.5,0.01719356],[7,2.5,0.01495233],[7.5,2.5,0.01264051],[8,2.5,0.01051631],[8.5,2.5,0.008530149],[9,2.5,0.006660428],[9.5,2.5,0.005018294],[10,2.5,0.003767823],[10.5,2.5,0.002747339],[11,2.5,0.001898165],[11.5,2.5,0.001264519],[12,2.5,0.0008248388],[12.5,2.5,0.0005535765],[13,2.5,0.0003488683],[13.5,2.5,0.0002137795],[14,2.5,0.0001238781],[14.5,2.5,6.509617e-005],[15,2.5,3.39542e-005],[-5,3,0.0002451456],[-4.5,3,0.0003979953],[-4,3,0.0006155684],[-3.5,3,0.0009329339],[-3,3,0.001345857],[-2.5,3,0.001988185],[-2,3,0.002849144],[-1.5,3,0.003884822],[-1,3,0.005174107],[-0.5,3,0.006660471],[0,3,0.008364709],[0.5,3,0.01031511],[1,3,0.01235795],[1.5,3,0.01450223],[2,3,0.0165354],[2.5,3,0.01844088],[3,3,0.01983436],[3.5,3,0.02082883],[4,3,0.02133358],[4.5,3,0.02128823],[5,3,0.02070817],[5.5,3,0.01979863],[6,3,0.01821609],[6.5,3,0.01626259],[7,3,0.01424764],[7.5,3,0.01213086],[8,3,0.01019035],[8.5,3,0.008268859],[9,3,0.006466097],[9.5,3,0.004957418],[10,3,0.003724719],[10.5,3,0.002745169],[11,3,0.001940411],[11.5,3,0.001308244],[12,3,0.0008651852],[12.5,3,0.0005758551],[13,3,0.0003675276],[13.5,3,0.0002385012],[14,3,0.0001492127],[14.5,3,7.314443e-005],[15,3,3.72384e-005],[-5,3.5,0.000191473],[-4.5,3.5,0.0003055675],[-4,3.5,0.0004799721],[-3.5,3.5,0.0007385293],[-3,3.5,0.00105315],[-2.5,3.5,0.001554785],[-2,3.5,0.002266804],[-1.5,3.5,0.003137288],[-1,3.5,0.004213018],[-0.5,3.5,0.00551402],[0,3.5,0.006971536],[0.5,3.5,0.008633199],[1,3.5,0.01031888],[1.5,3.5,0.01222394],[2,3.5,0.01411774],[2.5,3.5,0.01564666],[3,3.5,0.01703367],[3.5,3.5,0.0180793],[4,3.5,0.01851347],[4.5,3.5,0.01859006],[5,3.5,0.01817553],[5.5,3.5,0.01741717],[6,3.5,0.01626178],[6.5,3.5,0.01468736],[7,3.5,0.01292671],[7.5,3.5,0.01101026],[8,3.5,0.00922345],[8.5,3.5,0.00756762],[9,3.5,0.005978955],[9.5,3.5,0.004642918],[10,3.5,0.003535577],[10.5,3.5,0.002570006],[11,3.5,0.001824409],[11.5,3.5,0.001259681],[12,3.5,0.0008608194],[12.5,3.5,0.0005828783],[13,3.5,0.0003787368],[13.5,3.5,0.0002434085],[14,3.5,0.0001505519],[14.5,3.5,7.745292e-005],[15,3.5,4.095705e-005],[-5,4,0.0001357248],[-4.5,4,0.0002255649],[-4,4,0.0003600576],[-3.5,4,0.0005506782],[-3,4,0.0008093445],[-2.5,4,0.001201118],[-2,4,0.001722471],[-1.5,4,0.002406857],[-1,4,0.003227973],[-0.5,4,0.004253789],[0,4,0.005446336],[0.5,4,0.00679932],[1,4,0.008197852],[1.5,4,0.009746222],[2,4,0.01123769],[2.5,4,0.01263365],[3,4,0.01386307],[3.5,4,0.01486972],[4,4,0.01537166],[4.5,4,0.01555655],[5,4,0.01531142],[5.5,4,0.01466002],[6,4,0.01369118],[6.5,4,0.01242437],[7,4,0.01104347],[7.5,4,0.00943411],[8,4,0.007913624],[8.5,4,0.006556359],[9,4,0.005241971],[9.5,4,0.004096225],[10,4,0.00312608],[10.5,4,0.002266139],[11,4,0.001634431],[11.5,4,0.001181684],[12,4,0.0008147332],[12.5,4,0.0005360313],[13,4,0.0003551107],[13.5,4,0.0002228951],[14,4,0.0001345286],[14.5,4,8.036088e-005],[15,4,4.527891e-005],[-5,4.5,9.175488e-005],[-4.5,4.5,0.0001702452],[-4,4.5,0.0002787684],[-3.5,4.5,0.0004011895],[-3,4.5,0.0006036633],[-2.5,4.5,0.0009055028],[-2,4.5,0.001253258],[-1.5,4.5,0.001756835],[-1,4.5,0.002347103],[-0.5,4.5,0.003127987],[0,4.5,0.004042109],[0.5,4.5,0.005106617],[1,4.5,0.006192939],[1.5,4.5,0.007418597],[2,4.5,0.008585056],[2.5,4.5,0.009734198],[3,4.5,0.01068319],[3.5,4.5,0.01160795],[4,4.5,0.01216311],[4.5,4.5,0.01230431],[5,4.5,0.01219629],[5.5,4.5,0.0116391],[6,4.5,0.01086461],[6.5,4.5,0.009935535],[7,4.5,0.008893002],[7.5,4.5,0.007675618],[8,4.5,0.006469627],[8.5,4.5,0.005335872],[9,4.5,0.004239196],[9.5,4.5,0.003336417],[10,4.5,0.002581238],[10.5,4.5,0.00194901],[11,4.5,0.001429508],[11.5,4.5,0.001017616],[12,4.5,0.0006895923],[12.5,4.5,0.0004562329],[13,4.5,0.0003049306],[13.5,4.5,0.0001927053],[14,4.5,0.0001205748],[14.5,4.5,8.014286e-005],[15,4.5,4.763472e-005],[-5,5,6.032499e-005],[-4.5,5,0.000114748],[-4,5,0.0001958683],[-3.5,5,0.0002798341],[-3,5,0.0004227639],[-2.5,5,0.0006159834],[-2,5,0.0008340368],[-1.5,5,0.001175063],[-1,5,0.001604252],[-0.5,5,0.002135911],[0,5,0.002823031],[0.5,5,0.003610313],[1,5,0.004428428],[1.5,5,0.005360837],[2,5,0.006242996],[2.5,5,0.007136788],[3,5,0.007872204],[3.5,5,0.008579962],[4,5,0.009138712],[4.5,5,0.009257336],[5,5,0.009088407],[5.5,5,0.008715331],[6,5,0.008201539],[6.5,5,0.007474039],[7,5,0.006758425],[7.5,5,0.005975006],[8,5,0.005082452],[8.5,5,0.00420611],[9,5,0.003310479],[9.5,5,0.00262351],[10,5,0.002034788],[10.5,5,0.001534017],[11,5,0.001145809],[11.5,5,0.0008252413],[12,5,0.0005840532],[12.5,5,0.0003841919],[13,5,0.00025616],[13.5,5,0.0001605711],[14,5,0.0001009728],[14.5,5,6.846455e-005],[15,5,4.156762e-005],[-5,5.5,4.071459e-005],[-4.5,5.5,7.248987e-005],[-4,5.5,0.0001175548],[-3.5,5.5,0.0001754333],[-3,5.5,0.0002725037],[-2.5,5.5,0.0003915074],[-2,5.5,0.0005407106],[-1.5,5.5,0.0007630697],[-1,5.5,0.001073102],[-0.5,5.5,0.001467081],[0,5.5,0.001899732],[0.5,5.5,0.002422495],[1,5.5,0.003023389],[1.5,5.5,0.003678903],[2,5.5,0.004317995],[2.5,5.5,0.004929435],[3,5.5,0.005439891],[3.5,5.5,0.005974792],[4,5.5,0.006449226],[4.5,5.5,0.006568059],[5,5.5,0.006392343],[5.5,5.5,0.006178277],[6,5.5,0.005841099],[6.5,5.5,0.005366211],[7,5.5,0.00491371],[7.5,5.5,0.00440813],[8,5.5,0.003769309],[8.5,5.5,0.003128091],[9,5.5,0.002493047],[9.5,5.5,0.001977746],[10,5.5,0.001534838],[10.5,5.5,0.001155969],[11,5.5,0.0008337336],[11.5,5.5,0.0006096122],[12,5.5,0.0004446622],[12.5,5.5,0.000298758],[13,5.5,0.0001974886],[13.5,5.5,0.0001262647],[14,5.5,8.007965e-005],[14.5,5.5,5.097986e-005],[15,5.5,2.836026e-005],[-5,6,2.813023e-005],[-4.5,6,4.718639e-005],[-4,6,7.148542e-005],[-3.5,6,0.0001038899],[-3,6,0.0001559349],[-2.5,6,0.0002291988],[-2,6,0.0003405188],[-1.5,6,0.0004787663],[-1,6,0.0006712947],[-0.5,6,0.0009144391],[0,6,0.001201808],[0.5,6,0.001546483],[1,6,0.001954929],[1.5,6,0.002375743],[2,6,0.002777802],[2.5,6,0.003198453],[3,6,0.003560571],[3.5,6,0.003904834],[4,6,0.004191277],[4.5,6,0.004318669],[5,6,0.004246304],[5.5,6,0.004124766],[6,6,0.003921967],[6.5,6,0.003677133],[7,6,0.003448613],[7.5,6,0.003106221],[8,6,0.002672201],[8.5,6,0.002201206],[9,6,0.001761596],[9.5,6,0.001407846],[10,6,0.001043773],[10.5,6,0.0008115911],[11,6,0.0006150946],[11.5,6,0.0004439898],[12,6,0.0003138883],[12.5,6,0.0002098406],[13,6,0.0001371921],[13.5,6,9.075213e-005],[14,6,5.848986e-005],[14.5,6,3.420993e-005],[15,6,1.692253e-005],[-5,6.5,1.825595e-005],[-4.5,6.5,3.017732e-005],[-4,6.5,4.384929e-005],[-3.5,6.5,6.186057e-005],[-3,6.5,9.250914e-005],[-2.5,6.5,0.0001328797],[-2,6.5,0.0001933983],[-1.5,6.5,0.0002760926],[-1,6.5,0.0003815318],[-0.5,6.5,0.000523039],[0,6.5,0.0007330539],[0.5,6.5,0.0009513506],[1,6.5,0.001200232],[1.5,6.5,0.001447563],[2,6.5,0.001689532],[2.5,6.5,0.001970515],[3,6.5,0.002186827],[3.5,6.5,0.002399571],[4,6.5,0.002558368],[4.5,6.5,0.002626491],[5,6.5,0.002610845],[5.5,6.5,0.002588239],[6,6.5,0.00248109],[6.5,6.5,0.002366293],[7,6.5,0.002228867],[7.5,6.5,0.002038327],[8,6.5,0.00175832],[8.5,6.5,0.001442761],[9,6.5,0.00120348],[9.5,6.5,0.000960235],[10,6.5,0.0007295162],[10.5,6.5,0.0005418395],[11,6.5,0.0004117128],[11.5,6.5,0.0003110916],[12,6.5,0.0002186289],[12.5,6.5,0.000142793],[13,6.5,9.217743e-005],[13.5,6.5,6.268011e-005],[14,6.5,4.21289e-005],[14.5,6.5,2.212571e-005],[15,6.5,9.506848e-006],[-5,7,1.014732e-005],[-4.5,7,1.697972e-005],[-4,7,2.470876e-005],[-3.5,7,3.659359e-005],[-3,7,5.908607e-005],[-2.5,7,8.028452e-005],[-2,7,0.000106393],[-1.5,7,0.0001531847],[-1,7,0.0002287426],[-0.5,7,0.0002966488],[0,7,0.0004079498],[0.5,7,0.0005449046],[1,7,0.0006918609],[1.5,7,0.0008476598],[2,7,0.001028329],[2.5,7,0.00117028],[3,7,0.001276563],[3.5,7,0.001415556],[4,7,0.001490901],[4.5,7,0.001514479],[5,7,0.00153944],[5.5,7,0.001559063],[6,7,0.001504599],[6.5,7,0.001474975],[7,7,0.001400574],[7.5,7,0.001274775],[8,7,0.00110168],[8.5,7,0.0009093013],[9,7,0.0007424355],[9.5,7,0.0006016091],[10,7,0.0004594067],[10.5,7,0.0003443522],[11,7,0.0002631525],[11.5,7,0.0001992547],[12,7,0.0001419143],[12.5,7,9.417944e-005],[13,7,6.011386e-005],[13.5,7,4.031285e-005],[14,7,2.728503e-005],[14.5,7,1.36536e-005],[15,7,5.581337e-006],[-5,7.5,4.773432e-006],[-4.5,7.5,8.189544e-006],[-4,7.5,1.229094e-005],[-3.5,7.5,1.950959e-005],[-3,7.5,2.931298e-005],[-2.5,7.5,3.960803e-005],[-2,7.5,5.622156e-005],[-1.5,7.5,7.370241e-005],[-1,7.5,0.0001152187],[-0.5,7.5,0.0001619019],[0,7.5,0.0002339767],[0.5,7.5,0.0003138408],[1,7.5,0.0003890308],[1.5,7.5,0.0004658102],[2,7.5,0.000559575],[2.5,7.5,0.0006334299],[3,7.5,0.0007044766],[3.5,7.5,0.0007791389],[4,7.5,0.0008390944],[4.5,7.5,0.00084208],[5,7.5,0.0008594495],[5.5,7.5,0.0008934814],[6,7.5,0.0008675172],[6.5,7.5,0.0008488561],[7,7.5,0.0008328978],[7.5,7.5,0.0007694918],[8,7.5,0.0006560551],[8.5,7.5,0.0005456312],[9,7.5,0.0004577678],[9.5,7.5,0.0003789428],[10,7.5,0.0002777098],[10.5,7.5,0.0002154127],[11,7.5,0.00016913],[11.5,7.5,0.0001273782],[12,7.5,9.110388e-005],[12.5,7.5,6.14494e-005],[13,7.5,3.854077e-005],[13.5,7.5,2.350409e-005],[14,7.5,1.427351e-005],[14.5,7.5,7.570889e-006],[15,7.5,3.634888e-006],[-5,8,2.058569e-006],[-4.5,8,3.739088e-006],[-4,8,6.330346e-006],[-3.5,8,2.355903e-005],[-3,8,2.449824e-005],[-2.5,8,1.696285e-005],[-2,8,2.192415e-005],[-1.5,8,3.159442e-005],[-1,8,5.137874e-005],[-0.5,8,8.260823e-005],[0,8,0.0001290211],[0.5,8,0.0001831761],[1,8,0.0002126513],[1.5,8,0.0002403995],[2,8,0.0002782417],[2.5,8,0.0003205301],[3,8,0.0003611416],[3.5,8,0.0004002437],[4,8,0.000433357],[4.5,8,0.0004303379],[5,8,0.0004315114],[5.5,8,0.0004479643],[6,8,0.0004556575],[6.5,8,0.0004624506],[7,8,0.0004656567],[7.5,8,0.0004387536],[8,8,0.0003793194],[8.5,8,0.0003158119],[9,8,0.000273984],[9.5,8,0.000230187],[10,8,0.0001620605],[10.5,8,0.0001278126],[11,8,0.0001028206],[11.5,8,7.830204e-005],[12,8,5.660197e-005],[12.5,8,3.838953e-005],[13,8,2.380503e-005],[13.5,8,1.383482e-005],[14,8,7.978304e-006],[14.5,8,4.587591e-006],[15,8,2.524582e-006],[-5,8.5,9.657881e-007],[-4.5,8.5,1.95028e-006],[-4,8.5,3.504206e-006],[-3.5,8.5,1.115178e-005],[-3,8.5,1.237526e-005],[-2.5,8.5,9.83336e-006],[-2,8.5,1.13549e-005],[-1.5,8.5,1.580112e-005],[-1,8.5,2.551574e-005],[-0.5,8.5,3.92002e-005],[0,8.5,5.750337e-005],[0.5,8.5,8.228553e-005],[1,8.5,9.91564e-005],[1.5,8.5,0.0001138667],[2,8.5,0.0001307375],[2.5,8.5,0.0001503586],[3,8.5,0.0001754117],[3.5,8.5,0.0002063177],[4,8.5,0.0002176832],[4.5,8.5,0.0002035018],[5,8.5,0.000207976],[5.5,8.5,0.0002197901],[6,8.5,0.0002308824],[6.5,8.5,0.0002466464],[7,8.5,0.0002528225],[7.5,8.5,0.000237846],[8,8.5,0.0002083612],[8.5,8.5,0.0001782277],[9,8.5,0.0001418775],[9.5,8.5,0.0001100423],[10,8.5,8.653668e-005],[10.5,8.5,7.264077e-005],[11,8.5,5.856866e-005],[11.5,8.5,4.478147e-005],[12,8.5,3.252137e-005],[12.5,8.5,2.184672e-005],[13,8.5,1.339915e-005],[13.5,8.5,7.777622e-006],[14,8.5,4.573152e-006],[14.5,8.5,2.770881e-006],[15,8.5,1.615276e-006],[-5,9,5.515082e-007],[-4.5,9,1.219677e-006],[-4,9,2.156741e-006],[-3.5,9,3.363031e-006],[-3,9,4.709053e-006],[-2.5,9,5.690698e-006],[-2,9,6.253124e-006],[-1.5,9,8.109758e-006],[-1,9,1.236813e-005],[-0.5,9,1.797675e-005],[0,9,2.47704e-005],[0.5,9,3.346485e-005],[1,9,4.255791e-005],[1.5,9,5.08377e-005],[2,9,5.890541e-005],[2.5,9,6.777852e-005],[3,9,7.833792e-005],[3.5,9,8.979053e-005],[4,9,9.45669e-005],[4.5,9,9.313145e-005],[5,9,9.711491e-005],[5.5,9,0.0001040212],[6,9,0.0001109576],[6.5,9,0.0001174021],[7,9,0.0001215219],[7.5,9,0.0001182144],[8,9,0.0001043074],[8.5,9,9.201276e-005],[9,9,7.715406e-005],[9.5,9,5.848875e-005],[10,9,4.557908e-005],[10.5,9,3.918552e-005],[11,9,3.075648e-005],[11.5,9,2.30816e-005],[12,9,1.659559e-005],[12.5,9,1.098071e-005],[13,9,6.682769e-006],[13.5,9,3.901963e-006],[14,9,2.334602e-006],[14.5,9,1.45001e-006],[15,9,8.658412e-007],[-5,9.5,3.248671e-007],[-4.5,9.5,7.438508e-007],[-4,9.5,1.299193e-006],[-3.5,9.5,1.850239e-006],[-3,9.5,2.348328e-006],[-2.5,9.5,2.745881e-006],[-2,9.5,3.154707e-006],[-1.5,9.5,4.122928e-006],[-1,9.5,5.865323e-006],[-0.5,9.5,8.040028e-006],[0,9.5,1.072878e-005],[0.5,9.5,1.430755e-005],[1,9.5,1.857546e-005],[1.5,9.5,2.299166e-005],[2,9.5,2.742243e-005],[2.5,9.5,3.187871e-005],[3,9.5,3.619584e-005],[3.5,9.5,4.015978e-005],[4,9.5,4.278464e-005],[4.5,9.5,4.390961e-005],[5,9.5,4.533535e-005],[5.5,9.5,4.775974e-005],[6,9.5,5.040043e-005],[6.5,9.5,5.233904e-005],[7,9.5,5.269405e-005],[7.5,9.5,5.100618e-005],[8,9.5,4.718537e-005],[8.5,9.5,4.190069e-005],[9,9.5,3.534603e-005],[9.5,9.5,2.71142e-005],[10,9.5,2.127897e-005],[10.5,9.5,1.77393e-005],[11,9.5,1.37143e-005],[11.5,9.5,1.022677e-005],[12,9.5,7.299217e-006],[12.5,9.5,4.809277e-006],[13,9.5,2.932456e-006],[13.5,9.5,1.712689e-006],[14,9.5,1.017901e-006],[14.5,9.5,6.285792e-007],[15,9.5,3.756647e-007],[-5,10,1.59557e-007],[-4.5,10,3.683674e-007],[-4,10,6.345202e-007],[-3.5,10,8.576232e-007],[-3,10,1.015898e-006],[-2.5,10,1.21775e-006],[-2,10,1.606562e-006],[-1.5,10,2.253995e-006],[-1,10,3.133549e-006],[-0.5,10,4.269119e-006],[0,10,5.647953e-006],[0.5,10,7.071794e-006],[1,10,8.610117e-006],[1.5,10,1.058954e-005],[2,10,1.301801e-005],[2.5,10,1.550133e-005],[3,10,1.771732e-005],[3.5,10,1.971349e-005],[4,10,2.111592e-005],[4.5,10,2.148963e-005],[5,10,2.143981e-005],[5.5,10,2.149826e-005],[6,10,2.15542e-005],[6.5,10,2.161197e-005],[7,10,2.154276e-005],[7.5,10,2.125465e-005],[8,10,2.065819e-005],[8.5,10,1.906993e-005],[9,10,1.598837e-005],[9.5,10,1.210061e-005],[10,10,8.911344e-006],[10.5,10,6.755228e-006],[11,10,5.117238e-006],[11.5,10,3.8145e-006],[12,10,2.716026e-006],[12.5,10,1.802211e-006],[13,10,1.108003e-006],[13.5,10,6.421545e-007],[14,10,3.708128e-007],[14.5,10,2.218171e-007],[15,10,1.302674e-007],[-5,10.5,5.939211e-008],[-4.5,10.5,1.375441e-007],[-4,10.5,2.368123e-007],[-3.5,10.5,3.199517e-007],[-3,10.5,3.918598e-007],[-2.5,10.5,5.319167e-007],[-2,10.5,8.162968e-007],[-1.5,10.5,1.252473e-006],[-1,10.5,1.920157e-006],[-0.5,10.5,3.019602e-006],[0,10.5,4.260673e-006],[0.5,10.5,4.789419e-006],[1,10.5,4.645738e-006],[1.5,10.5,4.876109e-006],[2,10.5,5.893107e-006],[2.5,10.5,7.186595e-006],[3,10.5,8.285914e-006],[3.5,10.5,9.214796e-006],[4,10.5,9.905516e-006],[4.5,10.5,1.006634e-005],[5,10.5,9.797952e-006],[5.5,10.5,9.297572e-006],[6,10.5,8.604366e-006],[6.5,10.5,8.006683e-006],[7,10.5,7.795449e-006],[7.5,10.5,7.979307e-006],[8,10.5,8.314752e-006],[8.5,10.5,8.226698e-006],[9,10.5,7.191361e-006],[9.5,10.5,5.407573e-006],[10,10.5,3.653434e-006],[10.5,10.5,2.429904e-006],[11,10.5,1.684655e-006],[11.5,10.5,1.198666e-006],[12,10.5,8.391078e-007],[12.5,10.5,5.598685e-007],[13,10.5,3.468033e-007],[13.5,10.5,1.986883e-007],[14,10.5,1.10408e-007],[14.5,10.5,6.301073e-008],[15,10.5,3.589721e-008],[-5,11,1.626954e-008],[-4.5,11,3.784139e-008],[-4,11,6.6018e-008],[-3.5,11,9.33614e-008],[-3,11,1.290908e-007],[-2.5,11,2.078178e-007],[-2,11,3.601706e-007],[-1.5,11,6.170417e-007],[-1,11,1.145169e-006],[-0.5,11,2.184055e-006],[0,11,3.335304e-006],[0.5,11,3.554577e-006],[1,11,2.809516e-006],[1.5,11,2.262392e-006],[2,11,2.466164e-006],[2.5,11,3.00369e-006],[3,11,3.435902e-006],[3.5,11,3.733008e-006],[4,11,3.958295e-006],[4.5,11,4.030607e-006],[5,11,3.907656e-006],[5.5,11,3.580934e-006],[6,11,3.074849e-006],[6.5,11,2.619389e-006],[7,11,2.467472e-006],[7.5,11,2.65895e-006],[8,11,3.030965e-006],[8.5,11,3.249551e-006],[9,11,2.989995e-006],[9.5,11,2.270085e-006],[10,11,1.450295e-006],[10.5,11,8.450239e-007],[11,11,5.000279e-007],[11.5,11,3.171134e-007],[12,11,2.10882e-007],[12.5,11,1.394817e-007],[13,11,8.662414e-008],[13.5,11,4.909607e-008],[14,11,2.629811e-008],[14.5,11,1.426516e-008],[15,11,7.831502e-009],[-5,11.5,3.254276e-009],[-4.5,11.5,7.630869e-009],[-4,11.5,1.369958e-008],[-3.5,11.5,2.116599e-008],[-3,11.5,3.486278e-008],[-2.5,11.5,6.614921e-008],[-2,11.5,1.261958e-007],[-1.5,11.5,2.443782e-007],[-1,11.5,5.46815e-007],[-0.5,11.5,1.189942e-006],[0,11.5,1.903684e-006],[0.5,11.5,2.002463e-006],[1,11.5,1.44495e-006],[1.5,11.5,9.607321e-007],[2,11.5,9.170597e-007],[2.5,11.5,1.084515e-006],[3,11.5,1.204591e-006],[3.5,11.5,1.243813e-006],[4,11.5,1.26525e-006],[4.5,11.5,1.277832e-006],[5,11.5,1.245216e-006],[5.5,11.5,1.122563e-006],[6,11.5,9.139102e-007],[6.5,11.5,7.261132e-007],[7,11.5,6.764284e-007],[7.5,11.5,7.920044e-007],[8,11.5,1.006558e-006],[8.5,11.5,1.172341e-006],[9,11.5,1.13404e-006],[9.5,11.5,8.754691e-007],[10,11.5,5.413918e-007],[10.5,11.5,2.820651e-007],[11,11.5,1.377824e-007],[11.5,11.5,7.200651e-008],[12,11.5,4.275597e-008],[12.5,11.5,2.725857e-008],[13,11.5,1.683747e-008],[13.5,11.5,9.461957e-009],[14,11.5,4.921004e-009],[14.5,11.5,2.549141e-009],[15,11.5,1.348085e-009],[-5,12,4.74764e-010],[-4.5,12,1.130408e-009],[-4,12,2.137466e-009],[-3.5,12,3.785364e-009],[-3,12,7.598643e-009],[-2.5,12,1.645695e-008],[-2,12,3.368567e-008],[-1.5,12,7.209029e-008],[-1,12,1.824898e-007],[-0.5,12,4.248474e-007],[0,12,6.948243e-007],[0.5,12,7.304255e-007],[1,12,5.13449e-007],[1.5,12,3.179007e-007],[2,12,2.831286e-007],[2.5,12,3.25199e-007],[3,12,3.477023e-007],[3.5,12,3.363559e-007],[4,12,3.202634e-007],[4.5,12,3.140188e-007],[5,12,3.047227e-007],[5.5,12,2.710103e-007],[6,12,2.128583e-007],[6.5,12,1.633735e-007],[7,12,1.589919e-007],[7.5,12,2.113697e-007],[8,12,3.014538e-007],[8.5,12,3.773804e-007],[9,12,3.799712e-007],[9.5,12,2.980615e-007],[10,12,1.815911e-007],[10.5,12,8.801783e-008],[11,12,3.640058e-008],[11.5,12,1.485838e-008],[12,12,7.154605e-009],[12.5,12,4.158986e-009],[13,12,2.51441e-009],[13.5,12,1.401678e-009],[14,12,7.136911e-010],[14.5,12,3.564351e-010],[15,12,1.824956e-010],[-5,12.5,5.061484e-011],[-4.5,12.5,1.240092e-010],[-4,12.5,2.561736e-010],[-3.5,12.5,5.454587e-010],[-3,12.5,1.32032e-009],[-2.5,12.5,3.133892e-009],[-2,12.5,6.673219e-009],[-1.5,12.5,1.496287e-008],[-1,12.5,3.977527e-008],[-0.5,12.5,9.483966e-008],[0,12.5,1.564616e-007],[0.5,12.5,1.652857e-007],[1,12.5,1.174512e-007],[1.5,12.5,7.51618e-008],[2,12.5,6.924778e-008],[2.5,12.5,7.915975e-008],[3,12.5,8.183721e-008],[3.5,12.5,7.434197e-008],[4,12.5,6.542432e-008],[4.5,12.5,6.060696e-008],[5,12.5,5.713078e-008],[5.5,12.5,4.956183e-008],[6,12.5,3.78953e-008],[6.5,12.5,2.932162e-008],[7,12.5,3.204435e-008],[7.5,12.5,4.97383e-008],[8,12.5,7.830426e-008],[8.5,12.5,1.032128e-007],[9,12.5,1.06714e-007],[9.5,12.5,8.467899e-008],[10,12.5,5.134202e-008],[10.5,12.5,2.400681e-008],[11,12.5,8.956881e-009],[11.5,12.5,2.939062e-009],[12,12.5,1.054272e-009],[12.5,12.5,5.05023e-010],[13,12.5,2.878134e-010],[13.5,12.5,1.583606e-010],[14,12.5,7.949094e-011],[14.5,12.5,3.872393e-011],[15,12.5,1.937346e-011],[-5,13,3.963356e-012],[-4.5,13,1.023587e-011],[-4,13,2.430987e-011],[-3.5,13,6.400965e-011],[-3,13,1.798711e-010],[-2.5,13,4.509524e-010],[-2,13,9.676819e-010],[-1.5,13,2.13007e-009],[-1,13,5.534663e-009],[-0.5,13,1.307451e-008],[0,13,2.15861e-008],[0.5,13,2.311096e-008],[1,13,1.732426e-008],[1.5,13,1.279205e-008],[2,13,1.338282e-008],[2.5,13,1.557581e-008],[3,13,1.569281e-008],[3.5,13,1.352254e-008],[4,13,1.098771e-008],[4.5,13,9.38247e-009],[5,13,8.330693e-009],[5.5,13,6.924195e-009],[6,13,5.17951e-009],[6.5,13,4.240384e-009],[7,13,5.532268e-009],[7.5,13,9.937467e-009],[8,13,1.676941e-008],[8.5,13,2.281558e-008],[9,13,2.395821e-008],[9.5,13,1.914639e-008],[10,13,1.160017e-008],[10.5,13,5.343418e-009],[11,13,1.897689e-009],[11.5,13,5.447933e-010],[12,13,1.474408e-010],[12.5,13,5.153118e-011],[13,13,2.551202e-011],[13.5,13,1.360897e-011],[14,13,6.762243e-012],[14.5,13,3.253054e-012],[15,13,1.609715e-012],[-5,13.5,2.30342e-013],[-4.5,13.5,6.529597e-013],[-4,13.5,1.88248e-012],[-3.5,13.5,6.072555e-012],[-3,13.5,1.889188e-011],[-2.5,13.5,4.871282e-011],[-2,13.5,1.028721e-010],[-1.5,13.5,2.09902e-010],[-1,13.5,4.957986e-010],[-0.5,13.5,1.119515e-009],[0,13.5,1.84072e-009],[0.5,13.5,2.044551e-009],[1,13.5,1.755046e-009],[1.5,13.5,1.692653e-009],[2,13.5,2.091021e-009],[2.5,13.5,2.487593e-009],[3,13.5,2.455133e-009],[3.5,13.5,2.025635e-009],[4,13.5,1.526169e-009],[4.5,13.5,1.182029e-009],[5,13.5,9.629759e-010],[5.5,13.5,7.544537e-010],[6,13.5,5.556205e-010],[6.5,13.5,5.052242e-010],[7,13.5,8.04138e-010],[7.5,13.5,1.61398e-009],[8,13.5,2.842583e-009],[8.5,13.5,3.937922e-009],[9,13.5,4.170867e-009],[9.5,13.5,3.346788e-009],[10,13.5,2.028567e-009],[10.5,13.5,9.291936e-010],[11,13.5,3.233228e-010],[11.5,13.5,8.71463e-011],[12,13.5,1.964817e-011],[12.5,13.5,4.809219e-012],[13,13.5,1.808829e-012],[13.5,13.5,8.925797e-013],[14,13.5,4.383798e-013],[14.5,13.5,2.108244e-013],[15,13.5,1.046368e-013],[-5,14,1.013541e-014],[-4.5,14,3.340682e-014],[-4,14,1.2064e-013],[-3.5,14,4.565489e-013],[-3,14,1.509866e-012],[-2.5,14,3.942094e-012],[-2,14,8.116482e-012],[-1.5,14,1.489092e-011],[-1,14,2.983633e-011],[-0.5,14,6.14068e-011],[0,14,1.009455e-010],[0.5,14,1.246442e-010],[1,14,1.420251e-010],[1.5,14,1.909602e-010],[2,14,2.685317e-010],[2.5,14,3.231137e-010],[3,14,3.133216e-010],[3.5,14,2.492938e-010],[4,14,1.751883e-010],[4.5,14,1.22151e-010],[5,14,8.967145e-011],[5.5,14,6.54986e-011],[6,14,4.796917e-011],[6.5,14,5.019952e-011],[7,14,9.542501e-011],[7.5,14,2.060282e-010],[8,14,3.719861e-010],[8.5,14,5.2056e-010],[9,14,5.539841e-010],[9.5,14,4.455552e-010],[10,14,2.702147e-010],[10.5,14,1.235321e-010],[11,14,4.264583e-011],[11.5,14,1.120039e-011],[12,14,2.310756e-012],[12.5,14,4.323177e-013],[13,14,1.097764e-013],[13.5,14,4.53347e-014],[14,14,2.168748e-014],[14.5,14,1.05465e-014],[15,14,5.328708e-015],[-5,14.5,3.49186e-016],[-4.5,14.5,1.419874e-015],[-4,14.5,6.330915e-015],[-3.5,14.5,2.661892e-014],[-3,14.5,9.102899e-014],[-2.5,14.5,2.38888e-013],[-2,14.5,4.812969e-013],[-1.5,14.5,8.009783e-013],[-1,14.5,1.317984e-012],[-0.5,14.5,2.36639e-012],[0,14.5,4.069408e-012],[0.5,14.5,6.517696e-012],[1,14.5,1.100025e-011],[1.5,14.5,1.887314e-011],[2,14.5,2.82244e-011],[2.5,14.5,3.398557e-011],[3,14.5,3.248209e-011],[3.5,14.5,2.509547e-011],[4,14.5,1.658606e-011],[4.5,14.5,1.041237e-011],[5,14.5,6.801406e-012],[5.5,14.5,4.590458e-012],[6,14.5,3.37787e-012],[6.5,14.5,4.118361e-012],[7,14.5,8.967022e-012],[7.5,14.5,2.024845e-011],[8,14.5,3.708719e-011],[8.5,14.5,5.22005e-011],[9,14.5,5.570255e-011],[9.5,14.5,4.485968e-011],[10,14.5,2.721835e-011],[10.5,14.5,1.24356e-011],[11,14.5,4.280353e-012],[11.5,14.5,1.11297e-012],[12,14.5,2.213193e-013],[12.5,14.5,3.58334e-014],[13,14.5,6.220186e-015],[13.5,14.5,1.851099e-015],[14,14.5,8.245877e-016],[14.5,14.5,4.087833e-016],[15,14.5,2.133126e-016],[-5,15,9.850298e-018],[-4.5,15,5.074215e-017],[-4,15,2.650157e-016],[-3.5,15,1.18369e-015],[-3,15,4.118191e-015],[-2.5,15,1.083326e-014],[-2,15,2.156172e-014],[-1.5,15,3.382703e-014],[-1,15,4.819355e-014],[-0.5,15,7.921199e-014],[0,15,1.632077e-013],[0.5,15,3.748396e-013],[1,15,8.297885e-013],[1.5,15,1.571291e-012],[2,15,2.391722e-012],[2.5,15,2.869854e-012],[3,15,2.71391e-012],[3.5,15,2.050923e-012],[4,15,1.289397e-012],[4.5,15,7.347768e-013],[5,15,4.235302e-013],[5.5,15,2.607333e-013],[6,15,1.934308e-013],[6.5,15,2.723947e-013],[7,15,6.528505e-013],[7.5,15,1.515129e-012],[8,15,2.798858e-012],[8.5,15,3.952894e-012],[9,15,4.224813e-012],[9.5,15,3.405107e-012],[10,15,2.066694e-012],[10.5,15,9.440943e-013],[11,15,3.246085e-013],[11.5,15,8.40855e-014],[12,15,1.648469e-014],[12.5,15,2.506434e-015],[13,15,3.385395e-016],[13.5,15,6.531208e-017],[14,15,2.451912e-017],[14.5,15,1.236827e-017],[15,15,6.746042e-018]],"ignoreExtent":false,"flags":128},"119":{"id":119,"type":"light","vertices":[[0,0,1]],"colors":[[1,1,1,1],[1,1,1,1],[1,1,1,1]],"viewpoint":true,"finite":false},"120":{"id":120,"type":"background","material":{"lit":false,"back":"lines"},"colors":[[1,1,1,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"170":{"id":170,"type":"bboxdeco","material":{"front":"lines","back":"lines"},"vertices":[[-5,"NA","NA"],[0,"NA","NA"],[5,"NA","NA"],[10,"NA","NA"],[15,"NA","NA"],["NA",-5,"NA"],["NA",0,"NA"],["NA",5,"NA"],["NA",10,"NA"],["NA",15,"NA"],["NA","NA",0.005],["NA","NA",0.01],["NA","NA",0.015],["NA","NA",0.02],["NA","NA",0.025]],"colors":[[0,0,0,1]],"draw_front":true,"newIds":[197,198,199,200,201,202,203]},"115":{"id":115,"type":"subscene","par3d":{"antialias":8,"FOV":30,"ignoreExtent":false,"listeners":115,"mouseMode":{"left":"trackball","right":"zoom","middle":"fov","wheel":"pull"},"observer":[0,0,69.21198],"modelMatrix":[[0.4889189,0.6490813,59.09689,-6.484916],[-0.2521031,0.09550943,572.9715,-6.924098],[0.6033814,-0.4860439,191.5113,-72.37469],[0,0,0,1]],"projMatrix":[[1.915592,0,0,0],[0,3.732051,0,0],[0,0,-3.863703,-249.5012],[0,0,-1,0]],"skipRedraw":false,"userMatrix":[[0.5988006,0.7949587,0.09735678,0],[-0.3087619,0.1169746,0.9439187,0],[0.7389879,-0.5952795,0.3154975,0],[0,0,0,1]],"scale":[0.816497,0.816497,607.0136],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":1,"bbox":[-5,15,-5,15,6.746042e-018,0.0269021],"windowRect":[0,23,1280,680],"family":"sans","font":1,"cex":1,"useFreeType":false,"fontname":"TT Arial","maxClipPlanes":8,"glVersion":4.3},"embeddings":{"viewport":"replace","projection":"replace","model":"replace"},"objects":[120,170,169,171,172,173,174,175,119,197,198,199,200,201,202,203],"subscenes":[],"flags":5288},"197":{"id":197,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-5,-5.3,-0.0004035314],[15,-5.3,-0.0004035314],[-5,-5.3,-0.0004035314],[-5,-5.815,-0.00109626],[0,-5.3,-0.0004035314],[0,-5.815,-0.00109626],[5,-5.3,-0.0004035314],[5,-5.815,-0.00109626],[10,-5.3,-0.0004035314],[10,-5.815,-0.00109626],[15,-5.3,-0.0004035314],[15,-5.815,-0.00109626]],"colors":[[0,0,0,1]],"centers":[[5,-5.3,-0.0004035314],[-5,-5.5575,-0.0007498959],[0,-5.5575,-0.0007498959],[5,-5.5575,-0.0007498959],[10,-5.5575,-0.0007498959],[15,-5.5575,-0.0007498959]],"ignoreExtent":true,"origId":170,"flags":128},"198":{"id":198,"type":"text","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-5,-6.845,-0.002481718],[0,-6.845,-0.002481718],[5,-6.845,-0.002481718],[10,-6.845,-0.002481718],[15,-6.845,-0.002481718]],"colors":[[0,0,0,1]],"texts":[["-5"],["0"],["5"],["10"],["15"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-5,-6.845,-0.002481718],[0,-6.845,-0.002481718],[5,-6.845,-0.002481718],[10,-6.845,-0.002481718],[15,-6.845,-0.002481718]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":170,"flags":4136},"199":{"id":199,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-5.3,-5,-0.0004035314],[-5.3,15,-0.0004035314],[-5.3,-5,-0.0004035314],[-5.815,-5,-0.00109626],[-5.3,0,-0.0004035314],[-5.815,0,-0.00109626],[-5.3,5,-0.0004035314],[-5.815,5,-0.00109626],[-5.3,10,-0.0004035314],[-5.815,10,-0.00109626],[-5.3,15,-0.0004035314],[-5.815,15,-0.00109626]],"colors":[[0,0,0,1]],"centers":[[-5.3,5,-0.0004035314],[-5.5575,-5,-0.0007498959],[-5.5575,0,-0.0007498959],[-5.5575,5,-0.0007498959],[-5.5575,10,-0.0007498959],[-5.5575,15,-0.0007498959]],"ignoreExtent":true,"origId":170,"flags":128},"200":{"id":200,"type":"text","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-6.845,-5,-0.002481718],[-6.845,0,-0.002481718],[-6.845,5,-0.002481718],[-6.845,10,-0.002481718],[-6.845,15,-0.002481718]],"colors":[[0,0,0,1]],"texts":[["-5"],["0"],["5"],["10"],["15"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-6.845,-5,-0.002481718],[-6.845,0,-0.002481718],[-6.845,5,-0.002481718],[-6.845,10,-0.002481718],[-6.845,15,-0.002481718]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":170,"flags":4136},"201":{"id":201,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-5.3,-5.3,0.005],[-5.3,-5.3,0.025],[-5.3,-5.3,0.005],[-5.815,-5.815,0.005],[-5.3,-5.3,0.01],[-5.815,-5.815,0.01],[-5.3,-5.3,0.015],[-5.815,-5.815,0.015],[-5.3,-5.3,0.02],[-5.815,-5.815,0.02],[-5.3,-5.3,0.025],[-5.815,-5.815,0.025]],"colors":[[0,0,0,1]],"centers":[[-5.3,-5.3,0.015],[-5.5575,-5.5575,0.005],[-5.5575,-5.5575,0.01],[-5.5575,-5.5575,0.015],[-5.5575,-5.5575,0.02],[-5.5575,-5.5575,0.025]],"ignoreExtent":true,"origId":170,"flags":128},"202":{"id":202,"type":"text","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-6.845,-6.845,0.005],[-6.845,-6.845,0.01],[-6.845,-6.845,0.015],[-6.845,-6.845,0.02],[-6.845,-6.845,0.025]],"colors":[[0,0,0,1]],"texts":[["0.005"],["0.01"],["0.015"],["0.02"],["0.025"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-6.845,-6.845,0.005],[-6.845,-6.845,0.01],[-6.845,-6.845,0.015],[-6.845,-6.845,0.02],[-6.845,-6.845,0.025]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":170,"flags":4136},"203":{"id":203,"type":"lines","material":{"lit":false,"front":"lines","back":"lines"},"vertices":[[-5.3,-5.3,-0.0004035314],[-5.3,15.3,-0.0004035314],[-5.3,-5.3,0.02730563],[-5.3,15.3,0.02730563],[-5.3,-5.3,-0.0004035314],[-5.3,-5.3,0.02730563],[-5.3,15.3,-0.0004035314],[-5.3,15.3,0.02730563],[-5.3,-5.3,-0.0004035314],[15.3,-5.3,-0.0004035314],[-5.3,-5.3,0.02730563],[15.3,-5.3,0.02730563],[-5.3,15.3,-0.0004035314],[15.3,15.3,-0.0004035314],[-5.3,15.3,0.02730563],[15.3,15.3,0.02730563],[15.3,-5.3,-0.0004035314],[15.3,15.3,-0.0004035314],[15.3,-5.3,0.02730563],[15.3,15.3,0.02730563],[15.3,-5.3,-0.0004035314],[15.3,-5.3,0.02730563],[15.3,15.3,-0.0004035314],[15.3,15.3,0.02730563]],"colors":[[0,0,0,1]],"centers":[[-5.3,5,-0.0004035314],[-5.3,5,0.02730563],[-5.3,-5.3,0.01345105],[-5.3,15.3,0.01345105],[5,-5.3,-0.0004035314],[5,-5.3,0.02730563],[5,15.3,-0.0004035314],[5,15.3,0.02730563],[15.3,5,-0.0004035314],[15.3,5,0.02730563],[15.3,-5.3,0.01345105],[15.3,15.3,0.01345105]],"ignoreExtent":true,"origId":170,"flags":128}},"snapshot":"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABQAAAAKRCAIAAACEG5PIAAAAHXRFWHRTb2Z0d2FyZQBSL1JHTCBwYWNrYWdlL2xpYnBuZ7GveO8AACAASURBVHic7N2LYesqs4ZhSqEUSqEUOqEUSqEUrRkGIVm+xCuRr3qfs8/eiS3bSuLf8ReGGTcBAAAAAHAA7tUnAAAAAADAMxCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAAAAAcAgEYAAAAADAIRCAAQAAAACHQAAGAAAAABwCARgAgONyZ159RgAAPBC/5wAAOC4SLwDgUPi1BwDAcRGAAQCHwq89AAAOivQLADgafvMBAHBQ7P4FABwNv+0AADioTeglAwMAvh6/6gAAQEcGBgB8N37PAQCAjgAMAPhu/J4DAOAQ7tnuSwAGAHw3fs8BAHBQ53GXAAwA+G78ngMA4LjWiZf0CwD4evyqAwDg0BiDBAA4Dn7bAQCAFyB4AwCej986AADg2Si9BgC8BL9yAADAU9F8CwDwKvy+AQAAT0UABgC8Cr9vAADAUxGAAQCvwu8bAADwVARgAMCr8PsGAAA8FQEYAPAq/L4BAABPRQAGALwKv28AAMBTEYABAK/C7xsAAPBUBGAAwKvw+wYAgMPJOccYX3gC68RL+gUAPA2/cgAAOJyXB+Cp5V7z2tMAABwKv3UAADic1Lz6LAAAeDYCMAAAhxNjzDm/+iwAAHg2AjAAAIcTQiilvPos3oI78+ozAgA8EK/yAAAcjve+1vrqs3gLJF4AOBRe9AEAOBxS38C3AgAOhRd9AAAOh9Rn+D4AwNHwug8AwLGUUrz3rz6Lt8DuXwA4Gl7rAQA4lncYAvwmNqGXDAwAX48XegAAjoUAfAMZGAC+G6/yAAAcS2pefRZvigAMAN+NV3kAAI4lxphzfu05vO0QJgIwAHw3XuUBADiWEMJrA3ApRc7hhScwnMddAjAAfDde5QEAOBbv/WsXYN9qE/I68ZJ+AeDr8UIPAMCxSMx7bQBOKb1PAJ5Ww5BefSIAgIfjtR4AgM9WV8osz9KKfCrHvDzpvcMmZADAMRGAAQB4qvvzapyFFT8b65Z+ZRw2bjvuza614+VTedBXfQcIwACAVyEAAwBwy+55dR1Z10eeR1ZbszXjocfJ/O8XIg9hN7ctuC9Mwi/fhPzOLq7PU6QNAHvhlRQA8FXuz6sjsr5nXt3dJnba9+clSfjlm5Df1sWUS5suANgRL6MAgBc7z6sjsv4xr94uCX7zvLq7G9lJvl75btj3U75Fj07CpLiL7Nuy+eYwqAkA9sVrKADg/9zIqz+WBJNXX+ie4PSEJCz3KXe++91+uvHTIQADwEPxGgoAX+4vW1jJq19DvvnyY7r/+E0S3rFn1VsNAX5DBGAAeCheQwHgvZBX8Qi/jp3WNGvHJEwAvo0ADAAPxWsoAPzJW+XVEVlf/V3B2/l77FwnYfn3r5OwPY3/cibfjQAMAA/FaygA3MU2LpJX8aHsSbjLXZ0PUvqv5/C+BdXfhwAMAA/FaygA3MWSrbz1DG2Y6qtPB/g/j4id50n4nv9pyJHv87+gN/wLFAEYAB6K11AAuItvM1StM5Ct91LGiQ/y0D/c/NdIYX86jvi13rAemwAMAA/FaygA3GXzrt3e7v+i/hN4iafFzvUfiS4OUnqr/PaG9di3E+9bffcA4BPxMgoAd7mYH8Z7/QdNTAX2IsHpyX+puThSWC58qyHAb7UcbS5G3NF04PnnAwBfhldSALjLjbee6+6471ZOCZgXZqd1Era99K86k3PP/7sAAOC1CMAAcJd78gN10XhP1sP81Wex/0jhv2NNFQCOhtd9ALjL/W+Ubb1LjqcuGm/i70OA97XXSOE/kv95vtVyNADgCQjAAPCzX2xcpC4a7+PdAvAwBilZEn5y6cTbflve38VOXRsvOTEA+BEvTwDws7907pHb2vt7FoTxKql59Vnc8ruRwn8kj0IA/oWL+ZbEC+BT8GoFAD/7e6nkqIt+YcEnDutNNtze479GCv/RB31b3ocFXQIwgM/FqxUA/GyvvYLUReMlPvHPLpaEb48U/qNP/La81ki5m7hL+gXwQXjBAoCfWWrd8Q6pi8YzveG02/9ycaTw3336t+WFzgMwu38BfApepADgZw9qlmNv622Bi5UoPM7XJL1NEv7j/2qIar92ewWYbyyAd8YrFAD87KHdYqmLxqN9XyD5+0hhuYfv+7Y8zY/fOr63AN4WL08A8LPnjEuhLhqP8Jce5u/v10l4930Nh0IABvC5eHkCgJ89c4oMddHY10GS3v+OFH7nIcCllDcvWScAA/hcvDwBwM+eP0aVumjs5Z2T3iPcOVJYDnjb/2W9/57tHzf9EoABvC1engDgZy+cF0pdNP7o+X++eRO3k/A7DwGWE/6sALy5hPQL4J3xCgUAP3v5e2XqovFrL3/2vpyNFLY/JI2/Jb3zKuv7B8iLZ8gYJAAfgRcpAPjZm0SITV302759x1t5k2fvmxiDlCSnvWdVhZzSEfZsA8CrEIAB4GfyfvSt3iiPBeH3fAePt/LOS50vtO9I4R0dbc82ADwZARgAfvZuAdhQF417EIBv+/tI4X29c3cuAPgCBGAA+Nk7RwjqonEbezLv9CZJ+B1COAB8MX4pAsDP3jkAD9RF45w8KwjA/2szUlg+fub//D/i1QYAPhe/FAHgZx/0lpS6aKzRUekv7hwpvK/3n4EEAB+NAAwAP/u4NTTqomHoqLSLZybhj3u1AYDPwossAPzsc9+SUhd9cPLTJwDvaDNSePckbH+32vEOAQAbn/qWDgCe6XMDsFnXRdNg9lBoKfw4m0FKuyRhVuwB4NE++y0dADyBvM2VN7ivPot9UBd9NLQUfoIdkzB/sACARyMAA8APvikAG+qiq/zfMchP+Zg/4pfYJOFf/OmBP1gAwKMRgAHgB5IfviwAm02jrFefzpNI9E1T8pMP0yF2Wn5QA/Nv8uuRwvy8AODRCMAA8IOvHyQzmvp8fV20RF83Ofm3fCwB2D74bp++ff3TrZOwjRS+fTw/LwB4NF5nAeAHXx+AjVVvyvvvr6yLLlOxVd9R/CwfyCVy+StP68HkZ0qgehP3DFL6oJ+XnPxxykYAfJnPeJ0FgBc6VF/Wb62Lluibp7y5UC6RDPyK03mSg/zt5rPcSMIf9FJDsy4An4sADAA/+KB3pTuSt+lfUxddpuKu/L777kLoYz51P8UYKTyS8AcNbaZZF4DPRQAGgB8cOUWMuuh7ti++rYvLv+a7C6GP/NT9LI8YKfxQNOsC8LkIwADwA1uZefVZvNJH10XfWP41aUrf2hGaOtWP8ylJ+FP2KgPAOV6/AOAHBOBh1EW/81vzjRvLv8b6Yz3rdJ6KOtXP9feRwo/zrZPhABwEARgAfvBu7z5f7oPqou8Jt3WqbnKjO/Q3kR/Qp/ydAtf8eqTw41BaD+CjEYAB4Adv8qbz3bxzXbSkWflH0m+YQpx+fqduE5K+rxsWGzW/yf+OFH4cimIAfDQCMAD8gAB827vVRVvulUzrJif/llh7OwPLtXKk/Pv7MjAbNb/S+SClJ/+Zg5dEAB+NX40A8APqSO9hddG+ee2bY8mx1tXZdv/e7vNsS8SWme3Ib8rABODvdmOk8ENRWQDgo/GrEQB+QAC+38vrokdL5/XuX4nB54vAEnfHkq9FX/tAjvyODEynouM4Hyn80Jcs/rAC4KPxEgYAP5A3lATg//WSuuh1O6t18+eLi8Cb7tBj3diC8RdMBqZT0QHVZtRiPOJ/fXL/BGAAH42XMAD4AfV+v/bkuuh1pt3M/o1TXC8Cn68JpymNhd/1x5+LAHxwDxopzPMKwKcjAAPADwjAf7Spi37QN1MyrRU/T5fyrS0Cj1lH5wvCtg14fPwFk4Fp1Quz70hhuSsCMICPRgAGgB845wjAu3hoXfQ6027Km81YBL62JXiEXiul/vQqaFr1YmOXkcI8rwB8OgIwAPyADW/7ekRd9KbC2V367TYi7rWm0Outv18wEum1o2Lxzv4yUpiKGACfjnd1APADAvAj7FgXLal1nXgvLvAaG3o0Sp3Prx2hd10R/aEIKvjRGKRkSfie/yVSEQPg0/GuDgBuoeXpo40F4V/XRUum3fRzPq9/NpJvb9Q2r3cRryuiPxRBBfe7c6Qwr4cAvgCvYgBwi7zhY5jqE/ylLnqMPhqfXjtSIu6Nazeh99OHIRFU8Au3RwrTAhrAF+C3IwDcIm/+CMBP84u66E3B843656ll2hvrupsAnKZ0467eHH+4wR9ZEt6MFCYAA/gCBGAAuEXe80kee/VZHM79ddGbgucb9c+2VfhaByyzHpX00VXQ9neEV58FvsQYpOSck3/v3sIdAJ6JAAwAtxCAX2hdF31tpO2mpPlGhbNtFb5d2Ly59nOroFmpw45qMxaErTqaHuMAPhQBGABuIUi8g2t10ffXP9to36kVNt+Yb7RZQL598DtLzavPAp9h5Fv5H5o9c+RFLzSusT9CyadyuR1gH7/6xAHgNwjAAHALAfh9nNdFb+qZb9Q/j0HBt+cbbRLv5w5Dkm8RC3SYboZbW8u9mG/l4NJc24dvRz73SwGAfRCAAeAWAvC7GdsR9e17urf+eWzuvb2zd7OG/LnbgAnAB/Ffi7ebfPuXKVk8wQB8LgIwANxCKenbCkXe1C910TfqnzfTj27s7N0s+Y7C6Y8j3xaGAH+683Brf/r54+LtLniCAfhcH/l7HQCehgD8tmxRd9RFu+hSufyTGvXP5kal9PmS7+2u0W+LfPL+7l+8Pc+3L//hMmUawOfi9QsAbmGr23va1irX6rIbjbI2B6+HG003W1vZkm/O08gXH9oImnzyWhcXby3fvnzx9u+Yjg7go/ELEgBuYavbe9qs4o48LMlBfmTrftGb+ufpZmsruYWk5Zjq2Pf9iY2g5QsnnzzUj52l3nnx9u/ojADgoxGAAeAWAvB72mTaTR62umjJIRpLStjsDb7W2qqUyeticci1hDDZj13S77WtxW+L4dV/8d2Lt7uwb8irzwIAfokADAC3yHtcAvC7Oe93dbFVlXbGytkFd14XfV7YLNdL+pXwYllal4Lbp584CYkFuhvuHAt0cfHWhm+BPwsC+GgEYAC4Rd4E86733WwaU/3Y/3nURY8BwpsVY0mLY8V01DzLf+TyT5yEdNgFunW43eRbFm93RIs1AB+NAAwAt8hbPQLwu9ms90r6vdbVed3/edRF69JeDmNnb87Tul5Y7sqWfHU1OGhd9MdNQvrWBbo7F2+v5dtXn/73oMUagI/GSxgA3MJax7sZ673yb21Y1f69bvK8dn6V1UV7TUzWKEtLndfhSNKvJF65oSRku/bjGkF/YgAe4fY837J4+1a04zoBGMAn4yUMAG4hAL8bW++VUCofSLiVmGp59Tyjnvd/7pfnKaTiisSq2FJVHMuDdrdyb/YQ8oELcuTVucHv6Q3LFli8/RrsMAfw6QjAAHCLvC8nAL8VybRBa5N71bIlVfn3+U7ddf3zkJLeOOXqqm9tn6tzycJXyHq3NgrYDpaPNSpLWqufNAnpyX+1+d3irWDx9hMddoc5gK9BAAaAWyj2eysWdNdtmccG4E1fq+lS/bMNOppastUg3bb4tj3A1Wc/6qLXNwyxulC8BOHPicD7Pmn/uHhLvv0yn1hgDwBrvLEDgFsIwG/F6pPXl2iT5xZWy1TWV12sf15v99UgHfWGfaNvje1j7RftvLN+0dYd2lUv/+i68Sdk4P/donlt8VbcWLy1fEu4PSB2hQD4dLyxA4BbCMDvw5Zt14u6mwFI625V5/XP6wQrb+Bd6QfLv13M66LOUENIEv+8/CNJT7Kxrgm3nPz+GVhyaZi7WrN4i92xKwTAp+ONHQBcJe/zvP+wGbBfzPozry/ZpNw0pfHppv7Ztv4uN4xLANba6VrWP2dtrKVJUPtF60Ko123CVZ8OeifvU/55Md9aoGXxFg/C3wQBfDpexQDgKgLw+7Dl382i7qb5sxxjVdCb+uex9Xd8Ku/hk4biNKb+xhamTShp/UCSFl3bHtx2P5bnPCN+vXgr/5ZLyLd4BFpAA/gCBGAAuGpdTYrXkkQqmXbT5up8l6+1wtqsDFuzq+Wu2tUWgEfrLKtwniwex7zus6X10jn62vpjafiUD05O4xdu5NvzcLspTr4dbu3gP54ecBEBGMAXIAADwFUE4PfhdNX25HfWZgOwsVZY65XhnKfNO3a5I4mQduT6Pi0n60Zf+blPJz/3UHTasNxEwqfkUNdWhLVI+lIWPQ+3tir7nJ23NOnF4/DsAvAFCMAAcBXLHW9CN+We1T+PAUhrOrn3dKvwuvPztMrD2vvqdEk5JT1Yrh2l1MtVNY8TsLXiloST1UVvipNHuD3Pt0+oTJYHLesvGNiPPLEJwAA+HQEYAK4iAL8JW9T9ccyvWbfCajN+T+9qzsMWgNdXycG2OHwegMfCsn06aqpt5tBbtU1mSg0eh2cXgC9AAAaAqwjA78BKnTdx9+KYXyMHjwLmTdNm+Xi0sLJNxesbysF2rTXcWl9laXlcuM7V79YU993OB9+EZxeAL8ALGQBcZcWrrz6Lo5M0ex5W18u8G2Op1ro9r8W45GG5w/VWYdv9O5Z2NwOHNRJXPyYnjY5Zb9gnnIiCB7GpYK8+CwD4K35NAsBVBOCXszFFm67O09zt+eLxI9mu464Z2dDuUO5kBGCLvvJv+4FvFpxtTdjlKMF7ffy7tUmT83m3QI6vIQHYtrjLqyL7zAF8LgIwAFzFRJmXs6C7mfc7XRqAZGxl2JLtZilUfpKjnt0WeEcbrVHSLMfYB+tsbLRldF6Kqy0tv1uR/LudD76JvR5K9JUPRhJmSzCAj0MABoCrmPnxWqNP1T0DkIwFZu0arZN1T65qrZv1g1E+baOA7aqSq3xWdL5Rv59NANbm0imNPlhWBf1ugfPdzgffJIQwXg9t0BdJGMAnIgADwFUE4NeSmCoB9TzuXhyAZCwqa3IuYf2GfN22aqwn2z1L9JW4rOXR8iY+93Xj84fQABzzZm6w9+9VI0DRPh7nYgtoucT+7CJJeJ2QAeBtEYAB4CoC8GvZRtzzLHptANKIyrls5xiNdtC2qdgulBgc5EDnapxDY63B6fbGsTg8aCOumH1dVoaTpub3eobwjMXj3O6vNpKwDcHmeQjgbRGAAeAqeRtHr5dXGUl105B5ur4BeERlL6m2nnR4Xre/GnG61uKrk8S7vpMUSsu+cwCel7ysrNqXpQ9Wu9v3eqNP8MCD3N/vbSRhqvEBvCcCMABcRQB+IUuq6wVbc2MDsK0MW7XzehNvjEv7q/XqcY1BcvKmZFi7W/mqj1JsoJK3yUgWgOWfcT4SjZ27UBT6QheLVIG/Y3s5gK9BAAaAq4gTL2RJ9Xze77UNwDapaJqrnUNJPidJr5JvR/urk/DcgvL58rJ2t3K1pBCK0+g77x7WAUsluFDWxdXyDCnl5Oav5ZzjGYtH2Lclvpv9+jC3cu3yex4CwAHxugAAVxGAX2Us/J7P+/1xA7C849X2VLHIbeWnJx+3/lZ6jMTaUubCZh3rW1w9u7eifbBy8MvKc4vU1ljL+ZPdxfL2+q16TvF2Hw+y4/by9bP0xjP2xmF3XgUAF/EyAQBXEYBfZSzznm/3vbYBuA9A6u2cdUHYkmqMfbpv9EUCrH5knyddH143tZofwGkJtA9L0G0jj3SFufr1TeS50Trf7vmF/4Wcj/f+5+OA/7fXi+F5QL22wHvtktv3QAAG8CNeJgDgKt5LvUrfzXu23ffGBmAJxilLIp3Grm0rb5ZLqt5X8U4bSusnkomdK6lstgqrtl1YAnL06aSPdIw1aiSWe3MlxKrLvqUU78P7RM6c851tioD/tdeL4d8D8I2DecUGcA9eKQDgKt5OvcSof269q+7aAGzBWH5c6zgqdxJz6Y17gnY003LlVvw8pVRcKPlsxlKLy9oHy5X1UnMtVY6XSzQzlyQZWB4oBAmc0QL2O6BNER6k/a1nn7/0PDoAs/sXwI94dQCAq3gL9RLLNKPJb+qTr20AtiG9unc3Lxda32a9pNU8WwFzL4mW/zrtbnUy71cObQFSArDXzcG9P5ZNUSo6BNjFVEPSPlhtd7EGYInEb9IpPDWvPgt8oR3/trJ7AGY/MID/xUsDAFzGjspXGSn3/g3AvgYtTj5djLW2VTVLnPW2RKzdrJxWQevKctCdvScBeJ541Pb89upoDcN2ca2huFSKD313saRN55LtKX4HO7YpAtZ2/NvKQwPw/14L4Jh4XQCAywjAL6HThtrvpv/bABx0Q+9mjSrX0jb56hKttcgqudqPtA9GCqGVSLebzcu/pje7Kmm9qVjvpBa5ynYXS+CUO7Cm0+8gyNkQgPEAD2oBfe2Svxz2vwcAOCBeFwDgsiKp6U2SzZGMfb/n232vbgAute3InTZv0UPUMmaLtX3pWCKrW0XllEqcxx3Ny7/9tmHyRRPuet0ravmzRmJbHG6TUbNG5ff4OwlNy/EgOz61dgzA94RbAjCAc7wuAMBlBOCXGPXP59t9r20A1uhb6vkbXS2K3iwme19ap+gelWstwWkAznmzjKuNorME3JNoq/XS7UkhIVoCsKQCeZLITderxC9EAMaD7Bgj9wrAf2+dBeCweF0AgMvoqft8o/55unsDsObPoIOJtvXPWdtc+eokqVr982RRdTrJqzbc6Hz5WOcJ67bHk2grATgWr4vDOm84WeC0rtLvEIB5r49H2LEFtLlzbO/vulsxExjAj3hpAIDLCMDPN+qf798ArDty21CiTQTVXb5O+1VpK6y5/nnSrYw2BakfVnOb93v2RlkDcMwxn/SV1ixdXIpVs/GUxtvrdsd/+LL3wJZ1PMgjXgmvTSo6D7rnh7lL7rlzADC8OgDAZQTg5xtFznduALblX5fjea16cDq8SJs9T2HUP/f+z3FV71xaPD5d/rVFXV1YLicdriwA51h81FXl8Q5bbv3y7EnFPh4kpcQrIYBvQgAGgMt42/dk6/rnOzcAawOqohOAN+2v5Oem445i1HXayWlytmG+c/3zkldrm/dby/rmNtlIA3BddbhKqQavNdUhaequbqy42n2/dhgSf6/BgzBeC8CXIQADwGU7jr7EPUb983TfBmBrXOVa1eSm95PVP8ulS6huU47GqKOlZDpGzbRl+UGPeBxSm6LkWh9p+beE36KbiqsPTm/jxoqrrRi/dhGYv9fgQWiuBuDLEIAB4DIC8JPpau5Uprs3AGslc87ubOFTEmx0vauz3HA95UhXhrN+1mYD26O2fcJxCa8jG+sUpamtAMvRcwDQCcBt6ddlNx64l0y/tBUWy3R4EDbTAvgyvKgBwGUkiicba7xju6/O2p20jPniBmBdmq3pav1zu1RuqCHWju67gPWzlFpAbh2xtEd0XHbz2rJunaqma20F7XLyvY+0nWRrr+OiXydvmwb8wj+Y8HTFI9RaCcAAvgwvagBwGYnimdZrvLYULJlTPrAlXF13Pd0AbMXMetVZeaaWIs9V0c6GAc/1z6NKudc5t45Ymq5jX70d45C0dVbu/bGq77OUdPbv5GvwwRXnw7pAwLYNv7ALlQ0lftnD40uxtxzA9yEAA8BlBOBnWq/xSuyUnJmmHi8l+upK7HSyumpJ1eV4vvNW65/jMkvJIuuk7+OXYUVWtGwrwmlq835D6Be29CsXlqyPOx/aq6k1BqcgAdj7mPI2AL9wGzAbNfEI7C0H8H0IwABwWQiBJbWnWQ9A6rt2Z5Y8rRDaLrFMqg2uYt68Oc/tUg2jtVr5tEbWcFL/bHQbcNAEqwFYjpVE2243xibVEF1tZzXvGLa16FKSZGznQlq1ztK5wGm1tfjpqFPFI/B3QADfh9+XAHAZAfhpRq/mMpU+tWhlLA63jKmZM/aG0dH5uvkRaR52bXhv1FJj5gAAIABJREFU624l8VVvXmLb7XtyZIrVSph7oyyvhc2hzvG7jTbqsdyS8bwWXWpOLsmVoS5BXa63o16yWsZGTTwIpfUAvg+/LwHgMmpKn2as7l7c7jsWh20L7jSv5crH29BXa3HBduLWont3tcK5tjZacRtNi4+2Zdf298qtUvTL4nNrGG37fte7e219OIe8CcCSEaxX9EuqoCWihBfuP8b3cs7xMgjgyxCAAeAyAvDTWLMr+Ucz7c0JwNabyrLedgBSW7NttdF5sg3AVWNtjk53AvvT4uSca+j7hy1XRy8hcpW9Wxstbf5cyzIauGV13Z+ckmvrw+sH79ODX1EFTaciPAiVBQC+D69rAHAZAfhpLOJqCrXS5pXNBGDJn5J7NStP2YWybE5sAVRC7Ch07oOLssRl52vYvo23Bd4WVrWwuXqf4rJ6m3sbLd17XPN6YdcCcCjet1lI4/7GIXO59FMRgPEIPK8AfCUCMABcxtLHc1jEHf9c2wBsLKxK2tQNwG7V1Kqt8CafR6FzT6cx1qS9rLZv422Bt7VutvvU48fq7bxerB25+likflUfy1Tk2RFsXNP6Lqe5FvrJJKWkF84gxpciAAP4Sry9A4DLCMDPYb2axyLwtQ3ARntZlTaPN6flnblc2oYYVedHTO3rxrqDURtr1Vym9b20G1vrZi2dtjFL9vlqXrA2kbZWz6cBOOTgXPTFl1Uj6L6eXLfdtp6AVr14BJ5XAL4Sb+8A4DIC8HN4jZu6B9haQG+u3Vyic3hbbyq3fmfeoqdG3PlH1teN56Cri8YlL/fS6p+nebU2xNo7b9nn87WTFTzbsF/LxnPHrJRTC8BRz+XkLE4+eBpa9eIR5HlFAAbwfXh7BwAX1Fr9S/r5HowNQLLey5vtvueX2OKqzSJyY4N2W/6d2tjedf1zH1+Uc180Tqu65Ll42jbuuhx1uNE0r96u/vBhW471o7myuQdgbYKVgvxfXp4kY/fvSwIw+9WxO55XAL4SARgALiAAP4cFYNvl++MG4HlBt01LGlte57ip9c/tzXqf6zv1XblyvQTgnLYNrkzrZOX6uCM7enWtnl6Oeq9zkyvJ1dqsK0bnk5xDlNOfw+68SLx88DRUK+AReF4B+Eq8tAHABQxWfQ4b/Gsf/7gBeIwy8hpntRWW5EyvHZynmpaNu72V9Kh/dm2bcXJLgfKqqlNulEpZArB+voRXrcouoSfctm68DsAhpFDcOH5eil4+eBqCCnaXtU06r4EAvhC/MgHgAgLwc0j6HUXOP24AHinPae7U0JdcqlmHIQVfU6x2bW/O3IKupWBdkc2+p9JVVtRq6jZLadkz7P06vFoA7nl5XmqWh/beB4nUoeiJzMeP/s9PbgTNcxWPQAtoAN+KAAwAF7D68QRW/2xrvD9uAB6Vy7lUybC+Bu3APK/6VuclAPtWBN1jcwu61tBK76oGPfi0/llnKcWc5gbQfRfvqvRdA3D1fYl3rmxuO5BdzvJwxVe3nMM8Cng1NvgZCCp4BFpAA/hWBGAAuIBQ8QS9Vnn++HwDsIbS2ahctl26GkaT6xflXHycbEiSr3GahwrP7a6sc5Xe0jLwTKcfJb24B2A7ejVf2EYE9yfCvLDrJ++UPpauNrtgM5bWA5BOZhQ/GM9VPEIIgQAM4CsRgAHgAkLFE9j0o/Hx+QbgdSS2SJlaSpW35aVmLT82MZY434+m0Wrbhcdybw/Amo+X33q6/NsaXOm2X6uR7i22libOFoB7KYBcqNOGtKWW0/HCvYFWjXLTJMfYXVnufWYj6NQ86cFwGLSABvCtCMAAcAEB+NGs/nl8ensDsEXZUvqUIm1FlaIk5PlQVyX02odJrijr+udpbt3ct/jOtJeVZOWWY7VGOs/ryadNnNuO3/ZRO1RnFpdgHcJ1hJJNGA4htQLq0abrtNPWY1GqikegsxqAb8WrGwBcwKrao63rn3/cAGxR1oqf+98lWg2yLhpL+HPL/YQavauWPsd6bA/AOhCp94LW+D23d+6NoMtctXzaw2oEYF3UdbqwG3IImoO9PJDeSdVc3UcKz62mnzkJKQT5QsqTHgzHQGc1AF+MAAwAFxCAH21d4Xy+AVjy4zoA+1aevATgtiLcZxdJ/Isn9xOcFj9vZhG1amU3gqm1v7KrdA9wnkcHT9seVrrj19e+eNxybwwuthVeeaCWopNVPNsC9bwr+XmNoClVxe4ogQHwxQjAAHABZaUPte7/PF3aACxxdERiC59W+dxrjFv+HJ2rxg9K7yen4uPJWvG4KkXbxzvZ+KX52tYHKy0F1etmVhaAQ9HHjW0LcYwpJcsGuqKco2bnOVfb407PnYREqSp2F9vz/NVnAQAPwW9NALiAAPxQkm+XwHllA/CIxBIpR2lxz3rtP9raqrjs4lj+1PuJsaY8tgovj1i91iq3cJuLtr8aVcMtYJf1+ax7WMmp6rzfYKXSmpVHNtATS6171pyrbYiSblF+4iQkAjB2xwsggC/Gb00AuID3f4+jw3Und3sD8FKQ3OuO2+W1xFTHLF/d2VtdCuXkfub2V5tUKFFZ66X1oxAlIM+30vPJ1QYaLUevtvBaAE4+91TcAnDM0Y7SAJxcb89VigVgW/t9ziSkUop/5tBhHAN19QC+GAEYAC6gsdDj6P7byd3YALyOxOuNtf3I0WpZsl91sc3g7dfmXoJ82u9ZP5fQa4+iNdLBrauja4haYr3+hbjawit361OOru0qbgu7QZ8cwY6Sh9f20VbxnLOtBLv5gic8g9iriUegrADAF+MFDgAuIAA/ju6qXf32Od8AvI7ENvh3ObKsZvnGqCXQtdhnWjUd/Oj/PGKyCkFn9U66qJtqrM6f/GzbcKN10fUoadb4LQEzp+DK6mAfaujXhtY+eu5SbblXPmu10s+YhEQAxu4oKwDw3QjAAHABFYAPYu2v1jXPtzcAW//nyVpTtS2+S28r+SGlaLFW7lZrmOdsvOr33D8fbaXlTso6HFtD6dAm+o4APPfBCtppOoUUljgQgnNuHYB9dX3j7zwMyeLzcyYh0a4cu+OvKgC+GwEYAC4gAD+I1T9bap0ubQC2HcLj09H+qh85ZvvmXFzwsdhuYY3F2Y365xhXfZjb59oxq3WWdjkmt8qmLbDqnuES+iZhE0KO2vtKk3PymwBsG4Y1JnvdhzyKnmsuY2iwLQI/GpvVsTv+qgLguxGAAeACAvCDWP3ziJq3NwCvp+n2Lb4jiUrycznl3rxK/l1C3ytstcdLH+ZWn2wBWDtaxVzSHI5th7EVWp8FYLl/ud6WjkfZdU3J+ZMArOvV3vW6Z7m2HWndsJ5QRhokqBOAsSv+qgLguxGAAeACesA8gtU/378BeL2NVkuU4/pzyaH6Nwq7B73Pue3y+NH1NlTtcm0Z3bK3Xjsm/c4F1drOahOAUyqhFVcnr32w5oXnkrMPyx7m9sU4zd5z6fMYg+TcMwIwf6nB7nhSAfhuvMMDcJm75KEP97g7/4V3O5/vYPXPtzcAryPx+CFoQF1FXN2k6KI2ptJmZU4yqq7Btug5z0hSuiLrl5lJevuxt9HC8ZxrrdnVei06h9yXjqPXpeO5pbNujwxhnLZ2z6peJyG5PgxpHKn7ih//JJInKlkF++LVD8B34zUOwGUX3wM9543RO7z9eodz+D4Sbm0jrn16vgF4WkViia8jysaSxxZfZeHSa8/nWIOvbQ1WC5yD9WM2ui1XQvL8uYZUX3v3q/XEXlsS1l2P6eSGraGzPPI2AOsqb0/peiJFQ/jYBiw3nI98xihgnqjYV63VtVbn8j8IOuED+Er84gRw2bU31k94w/3y9/TyFpApILuz+uf1Aq+k35E5xzEjEnu/RFlfYo7zOOA2GSm51D+bUl8cXlUgGy2QdtU+79XXfr5uPV/4LABr6XWqNtLXNg+PXK3dgZzbBODeB6ttA5YAbEdaFfRDEwRPVOzOWkDLU6s90zUJsx8YwJchAAO47EYK3Vx1sUDaPr121fnl4+PbRdfPycbkikdo223jZgPwybbb1QbgzfKpq/OIoxBsMtJ4ImhALfPicFuAXe4u57Ee26uv83ytxdPV+uwIwH2DsP1HkkDbPDxmGml/IOfkEe3MtX20RvNlG7C25ppD/Vg3fhBJJuEJnaZxJPK/gzEDqeqeAn2OyeuhXMiCMIDvQAAGcNmdAfjGx9fy7cW7uufg22e1I3mfRwDencbIeRivubEBeL2BVrf15rnNcqtzziGPamidAFz7Wm6NSf5Z7i7G5Hsc1Q5Y1edaljs9DcAh98XnZdJwC7S3A7CuRrcdwC0F64JvcssEGbmfh+ZT5rVidxdbQNuCsG+YkATg0xGAAVx2TwC+EVDvCbrXHu7HO3k0CcAsrO3LKpDXHZ5vbwB27mQDsK6v2q7ddt16o68EYg3ALcpGX046LztXctV6Zqt/Xre5soS9WtEKqY8U7rnYEvIcgEdvLQkApXXeGgHYp6wHzAE4u7juwvXQP6QwrxW7u9EC2haE5fnNgjCAj0YABnDZOwTgGx8/FAF4d1b/fOcG4LbJdwmnMcU2bqhd1FpjrTf6xqyzeO1udcfvuGXLrDYNWO425Lze5auHjVXdfj8agJfOWzGOTcJa4TxPDtZ4kJI86CjV9rG0JeK+DbiENJ47FqIfh3mt2N09fcVtQdhKo/kTDICPQwAGcNmdAfjaqKQbi7e39wDf+fFDUVm6O10jbcuwm0vWx4z1Ye9PcmMOfom1rSx5qY7WNtC9ILkv0o7ZwfMyce9TFWqrhm5v1nMfcbQuUE5Fc+wSrVcJ2QLwvPysZzICsPbSCpacne1PLqmMe7VOW49bJyMAY1/WAvq/jucZCODjEIABXPa7FeBrN/9x4fdG0B39tH444/0QgPdl0XdT83xtA7Ctmi5bfPOkM34tU7ZAuZ70K3co18YadHexBd4Ra+cnjO4aTpqZlxOwkFxO6qVzqa76fs+nCVmbdWlfNLtXrZCWuGtZWs/Itxrp5G3RuOYy7tUaaT1uhexGtSrwC7z0ATgCAjCAy/63VnlzCQEYQ+yNpe7aALwqPVZFbjQCZbub9QZg3X4bQluoTX3vrhU9r1JyTNUStQ00ao/Unkt25EwirgTgvlrb1pmnue7Zknlby+3d0UYA1ohrATgH63l1HoAfV01PAMa+1i2g93Ktpf/9h92YCwAAv8BLCYDLrr0RuXHJjyH2zstfHoDpLbSvHiAnt94AvAnAIxLbbt8R64oLS95tBcrra204sERVV8Lyvt1GJc238iVaou4BeOnyfDKnSFee7Xfi2Lk7J2S5lT6EruX2sUPyuHa2FoB1fbvEvrAcgt56vlvnHtgHizyAfe1eVH/nHpZ7fnfcvgcAuB8vJQAuu7G/99qRmwuvfXrPHuAb1z4BAXhHo/65r7428vHFDcCb+ueatea5x+aWLdf1z3K8VkfrGqgu3i4/MWuiNR8TarTMbP2c++puP48lXetWYfudOJaY+0Tgvl1ZbhdiCrEF4ORHgNfErpuE24imlpk1p5dq1/rV/uV9Ma0au9Mm5/s9WW//wfSew+68BwD4L7yOAPgAT37TQ2+hHdns3/MJwGM12NgqsQ0nWqKs/CCiGx9rx6vV8m3J0VZXLQAvtcB5qZru/Z9bBNUAXN1Jf61VI2hdqrUm1eslZi0ydi00e22jFaLXUUuSdsPI8xpxtXY69VrnEILG4WLXbr+i/dCrHLtzd7SA/q97+/GS+w/78SoAuBOvIwA+AAH4c2kHqVZ7vN4AbJ/a3l1tZGUF0rUXDC/vwF0rLZ4/nlonqrFAVaIfaVgD8EjUFjrbvcijh1jHGqwu1q63OM67fC1X66nmJWHLHZbgartIz1nvNaacrOfzCMDWZTrXeeZvzhqAQ4+8sYX0RwRVdqpjd/u+0hKAAbwnXkcAvLWXND4hAO/I6orXS74SeiVq2rJw7181udbeqnfA6tonPTbPpc/rAUhW/2xXSgSVmD0/pLPNvZa0raba1mDXg4jVvMvXUrfu9Y1987BlcvvUAnA7HX1iyFU6fHjqi8wWgEvN2pGr3acE4Ox7ALbBwo8oVX5EvyIcmTy3960p2D0Ak34B7IKXEgDYCrrXNL/6LL6BtbZa93zWblKTO98AbDXGFhc7Sb9xjs1ts+56A3AteaRhjaA5nsz4bXdk9c/a/9kWenP21S05eTxKqr14eQqpfdTTrxzZEnqa/+NcsB2SNcWRpTUYtxXgMdgpyhfikoXzlNYL0nvizzTY1+41BQRgAO+JlxIA2JIAvGMnmCMb6dcCsERKyY0XNwD7Glz1S1Oq1g5Lt+ya9sZ3PQCpxrDUP7vJMqp+sprxa/XPOknYFnpjDOUsAIcQfbG7bQODk0X0fljrvGV33s+oBdmaNe7aIrAG4Kxf5kjXbVE7WdW1Ze9Vt+ndEICxr4e2gL52yV8OA4Df4dUEALb2bYV6ZJZ1rYuyLatKEnZnv3r08jiFkiRJ9ovkRxD6qKExuGi9PVjrn9vPyFZ8L8z4dS7WMG6iFzsXqt8G4NRm+jax+JCzbVru17b8agG4lFWLoKQBWEumJQpLqk4pVv1vrmWysmeX1wF41WxrNzxLsS95Rr1tACb9AtgRLygAsKWtf3evWD2eZbRvi8G2rLouh14fpqu4udpSqi3/xjxXSsuPI5d5JlG/yfjEloX7iKN1kbQE15jHTXQNNvZ67PWj1zhXXecsqVbnCa9Pr60d2xlaAJ7PQGc7SaoN2esYpJR8ivKx10HBrexZvlbbS1z73KXdt+vyLMW+dn9G7RWASb8A9sVrCgBsES120Uf7tvQ4ukCfR1CtF869A1bMWiMtabU670KxIbrW0yosJc/TulezrfH2ALya65uTz27Z0hjbBuOlUnqmSXUeX5QkzebTnNpitzXTOpm727YTyxm66jWtS7Su2hpa7iEnbTo97rbOY5V274NFKsC+HvGMunOc+43DeJ4D2B0vKwCwte8wzMOyNlGWfpepuZs+zG3nrSRJm91rrZhz8BIgJYvqD6Gt8OrlqwBcQ+/VnJeV4NZuev5EU3fRCb69RrgUCae6NnsagHXBOMzZ1GlJszuvVPa+VA3ApUgg75lW7k13/NaqX1rOcv562jlIItYgnSbvlshrK8D7vo2X5yfBADvavQX04Gbnl/94mLvkEScJ4FB4HQGALd5j/Z1EUNvra12vbFet9ZfaHOna93tE2VSjdn+Oc6X03EDZRgT3fLqqfx6puE/xtcvl4rpqvhwlDkddHj4twO5ThdugpBL1PHvqXtOOaKnl2+yc3jaFkuTeWz8t/YpKsAAsR4Xa+nnVFoDnR7cOWOsJxn8nafxBcQXHxFhpAMfBmzwA2CIA/93o/GzNouzCixuArf557N4NXttQLZuB26X2X1sHDqnkstQ/j1SpDa5q7pdPPuVa3Nx8WUJpqXrw6JXV7rOv0bY5Rb7qnl6d6Dvf4Xy/oRY9JQvAVe436KlYAK6TfJB8C742G0mPz8V6bskH0xyA920ETVzBvmgqDuA4eJMHACcoLt2F9VK2deDRVPl8A7AeIOG09N27+m8XQ51HJc17esdKry6uasCsdvP1D2o0zbIH1Zv61HtktcVSXUAuSwBehipp0bJ2nLbl3G1MjdGSrbaLdmmEau353B4uVf0SLADr+nbrgqWrv87rCcwPtG8jaAIw9sXwcwDHwZs8ADhx0usIv2WlzmMd2JxvAI6lN2q2Tlc6kUj7P8+V0nPAXa/0+nmv4rrls/Vwtv299qB6U6s8npOuRuJSbXjvtNRW673LbeXELABvU0BKFoAlcDqXS1hmJlkA1jZXOrDYyZHy31o0b+t+ZufljCz32j87JlY5mbT7YCUcGJ3/ABwHARgATrReRwTgP9kMQBqXX9gAHLPlU/mWa+fkthtYJyBNbt3qeT0AyVVv5cTLEu6kiTDnPsFIw2ouelNrwTwnXY2gqQfgJTxnybTOPtFZTdVvc2U7tC0pBwnA4xGTLu+mcXoWgLUuuuoXo2eoTaKD7VvWScVzt+ldULCKfVH2AuA4eL0DgBO0F/o7G3qky7Cr3zKXNwC3+mfJh/L2u3j9j3bDylZy7Czglly1V7PdpLRhSdn27U7LkpVz1qt5ajF7Xd48UnS/lS1Nzwdo56rWunmaA/B2nbYl19apy4d1ANa14JMAXFPsddGalYsWS7ckrGOBw86TkIK25iq73R2OjRc9AIdCAAaAE7wX/Dvb96tBdJV4L2wAnjs/W52ydWPWC6tWF2su1YlDWhRdfbAdtJYwLUyu6591p24bBdyD98jGIYytt7YGawF46Q+dXa7JsqkF4O0Pvz2Y9uXyTgLuWB9eB2CdV1y1B3Uvww4hSlT2euZ9LHC7zx3X2ChYxY7YUg7gUAjAAHCC94J/ZBHUln/XiffCBuDYeztrHI2prQKXnlElP1ZXdCXYB8l6pSdh3WTb7kRuuCRVmxXcArD8E/Qn2C638bvzcRabrUF0b6mV4mazsXXPOtFuJvnWObduZGVr0faxnrN2z0o6v2mKNnZYs/DcEdrWfndsBE3BKnbElnIAh8JvUAA4QQD+I1vpbVtvf9gArJ2i6lykrJuAWwGzTQDWGl+XSrtkvp2O451nAXu/6ow1ZgW3B13qny0lz5XHlrTDFHwsdoBOB7YQO2dTubmfy60XTrfzulaePZ4auabRUNoCsJytjgSetNx5BGA5WyvwngjAeFdsKQdwKPwGBYATBOA/0iXWKf28AThPru1j1b5XIVv7K9sArGN+nbOEuc6ckldTjRZo2+HtqtUR2v5qipar2+e+L/vOuVNjc5WQWvSAnJeIPmfTVupcN+dZJRQn572fByqpUk4CsC+6AVgCcG+yFXJ07cRS6kvcdbdJSPRpw76oqAdwKARgADiRmlefxaeytV8Lorc3ALdmVqX3svLBFkl1yTZ7a91sJc3Lcq4VRU+6kTfHYrXTGgNXzaCt/rkvjo60ulp4taQaix6fSljOcM6mkml7PJ7pQzgXomt9p1YBeO65Zaftc5uB5HsALiFpU2vrfzWfwjo//wV/o8G+KCgAcCi85AHACQLwX0iO1G20U9gk3vMNwFqNnIpzS/uryeb92lCi2KcoLXE0zkvKKemgpHbfNm1oNIN2bYhwz4YjGK8WXjVv56h9qko5qdCep/RaAB6FylkDvS5Hx+Akc64DcK1ljBTW7cltCLDuFm4nqUe60uc7tce1DEwAxruhoADA0RCAAeAE/WD+QkuIW++rTeLdbADWUueoM5DkjbclzD7xyIWa4siuGhtT6TdJ837dUqrcst23rrDOFct943FKm43B69xpSVUC8Bga3M3HbAKwLf/qc8Lr7dajjDYB2CddoNY2Xe0L1yMlX899rmNcyrH/jr/RYEf8PQXA0RCAAeAE/WD+wtpQ2QykceH5BmDtC+Vrj4XOjz2yJbTG0XN2jatGUzqvqC0p6w+n9YuWj2tMI1LKQ2h8ze2B1luHV7lTA2lKsYZeTT3Mx+jW5fnnLxlTw3nbxxslWefTANwqtO1jzeGpfclOV7/lnpcA3FZ+U29xfTq7+Ld4imJH/D0FwNEQgAHgBOni12wN1oqfb28Abpt8NQpqpyibVNTWdV31NaelqVWcV1nbBuAxAEnLnu0texuSZAu2cpdp7MuVO1z/EOdtwFoCHbMVaW/P3pZqp+gl6ab5AjvLUiwAz0c1tY6Qr5XO7VSrd6F6i9Z6pLWhztkCsO0I/nsj6BACT1HshVc8AEdDAAaAE6SLX7POzxJTb28AtrpgK3PW5V/vrUQ4uLaIOm8AFq5VFcttc/IjbeqicS494s6zeS1yhzivym6a+sx5WIN3KLZGvT17W6qddICR3afOMYq96DrIA5UegPsS7ioAS6b1QR9ajtNW1VMf1KTNvdp9jYHEm2D+O/TsxY54OgE4GgIwAJxozX7Lq8/iI1nz56ktxo7EuymHnqz+2RZWY+5zekOoIfYJwPPUIqtitorilJzds40U6hXLLSrb7t0+fNi3ALyufzZt7dXqkSUAj9LlE3MADiVZb+ol7NbqnYu1t+kaYWEE+xGA257lYAG4N+hq7bX0AN9nGP+92tQ5R2LBXmgBDeBoeNUDgBPeewLwL1j9s4RV64C1vnyzAVhyoGTdaZrzodPmyX0C8JTH4q0tEWtgzD6UZQNwz7Z9gFK1LKzDh7NuKtawvRr827UAahf3kHyuZVP5f6/1zqvl36kt9s4BeF3DrAG4lslieVuslqtj7lm9B+C2uG0nGUKvif4jEgv2Uttz+9VnAQBPxaseAJygIPB3rAfV1CYh3dgA3Oufa9IyZmv+3NpDaXbdbABua62SqCVajjXkpYTYkvN8YSjJArPVIW9PrgVQPVyXkvPlANwWiOVUQ40ny7+Nc33b8DoA2+r0fPeT1W/LI4wje4PrOQBrs2j/1wAsT06G1mAvtIAGcEAEYAA4QQD+HUmVFnTHB+Py9QZgDbo6xFdXcmtMvZJ56r2pluG98zZeua011lpf2O9ozoEx1VFWraE0X3pDL+E2lpS97QG+cEBbJdaC7RK0gfMqFegqmV/F2vmL08eqeZoDsH6lSQOwBWz9UnwaK9WjC/Qf06sklrDLNGGgtYAmAAM4GgIwAJwgAP/O2Pe7jpfnG4DblKSgRzpXUy94trHAmnLnddel1Fny3nn98zT1PbVtNdbXYF21RCx+03HaaD+qoNOPbF7ROpPPR9QRgPs24/FVlOLCKtbOd68BuKTxdenic5KY3Y9s5dS573gupRdgu21/rv/Fkh12RAtoAAdEAAaAE+yI+4Wx0Xez43c7D6mFQf2/lLQ82DbF2uUlSNQdAXfJmTFKarW2UsuFVjjddu3a3uO+MKsDh/rBa5I85WKJphJZtXS6+gsBeNJQnYvmbe/qehOxLrr6Vayd716ru8tSsC05XPJwi89uGgHYUm8pVjs9rwf/7zd4wdRW7IjK+QzsAAAgAElEQVS/9wE4IN7nAcAJAvAvSLC0dVfJhOv8udkArNFwbvWsQdV7S7S2gbbGJQAvKdG5MbZ36cBso5Ja0XL7KNoOW92Cm/2m51Y/3GWJpnIyGkhLuDAGaWodwL0OHB6zhc2NACx51z7WL6EsAVhHN+U25dj6Vqe+RXmOw//5/T35Wliyw27oKA7ggHifBwAnCMC/MIqKtQx4lS3XG4Bt1VZbPeeQXdQJwHPMPR+A1H8IWu9rk4J9H4DU77dF4XZRKHLrYlN25dNSs6XlNQ3aOqdIb38rAMdYtMGU30zr1X2Svq/rrsuwrUO1fazjnErMNdV2ZNG66SnIf3Vfcpjm2cL2VfwlwBKAsRdaQAM4Jl74AGBBi91fsDVY+/jGBmDve/2zRMTeHrlFSRuA1DcAN2m0gta0pyXN2lhrtSt4icIhaF/ndgbWuUoedBOALYi2Hlt6K42gOV7cJ6w9q9qC82Zar55FjCMAjxZUmwAcSko1bgOw/KfdZtxnX6z+LWpWsRf2kwM4JgIwACwIwL+gM3hb2fONDcB9+Tfm1oPZ9W7IbSXTBiBp6+bVAKS+xum9rS3rInAsY1fwOLIWzZajv5TczA5en57WHjvdSFxDX2GWAHy+T1hnDjtt2nwegIN2kI52JraPd/kC8zK12Ge9XY/QU9KmWnKLUqxUewR4jcp/6OJMAMZeCMAAjokADACLos2KmDHzH6wHldU5b3b8rj+15k+hlTRr8bMNMWpBrk/QXZUdr+ufLUJLZnSh9PtdtgJPsQadUNSOD04XXW1s0voMdcU1to3B80PoUu1ZALYK6qmVQKdc188CzZxFxxFfDsDtZJIk+5TkfLQqu/olAM91z7o7uN0wpT8FYGpWsRfK6QEcE79HAWBBAP5f2l5q/lWyGfk7PrWNr5IsQ7U20E4nAJ8PQGqWUudW/2wBWPJkrKlfvUqAlpxz1DLj6LIt266nHFknqj5IaV7YPQ/A1ibaulRruXW5HIDLVGzkb7/VVHRB2yqcXZLcHIvWebexxG24sZv6UrIE4JDshm09+JffbSoUsCN5LpW/NGQDgM9EAAaABTWB/0vrmn/aANynH7nJNuK2VBvXG4AvD0DS3OmsW1XM887eVf2zlVjrxt1WVZx8ttSqu3HnHldybGpVze2c+uptLHnTKbov/7aGVZqWa1nHTFt0DdVvA3BtFd1tU69+IaHIMdXrtuSlbXXsG3/lcrvhSTev/8QfaLAjWkADOCYCMAAsCMD/xeqfrc75xgZga9HsYtY+zM6lWEcZs14uoTGMXb/zNfJp2wDcLwzzzt51/bP2vMpaZtz6a42V1XUAlkvGMvI0h9dUThplycV9+bed0CYAj065425HfC0pSES3Ds9tWFNxGoFDDFHrqNuYpxKS3a1+OfMN53bX//8N5/mJ/VBOD+CYeO0DgAUB47/o1tyxzHtlA7CVNGtfq1b2nLXiuYwQaX2hT6qa7UPdnRjWEVr/r5RN/XNs1cYaL2McK6vrh5Zku1REz0dsOkXrPYy5wylZvfSSckuxquMx67hfpS2fNQDLQ9jm4p7S5UxyklQvyVxH/vrYH8M5P6+3yfG/C8A6kInnJ/ag062pJgBwSARgAFik5tVn8THW9c+bDcDLZOC5RbPExR4Cw6r+ObY+UqOqeWwAdi7WsI7QugAb/ab+uS+ltkG70xws205f/VTrn6391dCOKPWkU3QvyTZFBybJnY+VZvubiI1ikn80rPsaqpydbhVuzaXjZOOFfbvbttNYF8aDl3zcA3BbI7Y2XdNyvv+NrkXYC3/sA3BYBGAAWBCA77epf15vAB71zzY1t8fatvyri6lzBbAui0r+jW7UP/cNwO0GI1H3e5ii5sp55dSWeftirNVYW+JNSwDWXcHFnYz8bQF4PSpJH24s/04aZJPWQ6exSCs5wcV2yaTBXW6r8TUEDcMt6NpDW3m1ftoCsFxbgtwu9y+5VXSPAPzrRtAEYOxFnkta708TLADHQwAGgAUB+H62HGoZdbMB2JpcTXOg1QQbtHtzDLXXG9s9eN0AvM6ffd211T+PRG3Xa6xMJ/XPKbd7swewzb0tWFqFc9Z8HN3m19zcCHpcriu063WwWkcAlnvQVl7yabtJakXP8oGu66ZiBc96P+2hlwDc0rucQCraGlr+yalK9NVJxq6n1xaHf/U9p20vdiLPpaAztBV/VQFwKARgAFiwwnYna/I8YuTYH2ssGeoO4RxTbeu02o2qVSzHXuXcG0f1Jsp6q2XCUat/HmvIlk9LjmPjruXtPtbX/mO1za0PlgVgjeDZbbo9j0bQ8sil1u3yb5O0r7RmXomw8oWEGOwpoRcWfSxJs64t9NpKsnbw6mfd6sCLnkxonalTKBKAay6SmeXMku+Dmqwm/BdJVhtj07YXe7DnUtX/FWRLwvztD8BBEIABYEEAvpON/z0fgCTh09KvxEWdppuzBsTaioW1nbK+79Y2U663Tc7RlVSsS5a3Hb5z/XMfgDRPRarejbrlbf3z1KufbRnWcml0WTforrYlt/PrraL1/uuyG3nNmm+1BVw7gf6U0NRdg61jjwJmTby+V1n7loZL0WVo3Z88Sd6tUa/MJVf7wi1iWDeuXwRg2vZiL5vnUtXahyQxWJ7wVBkA+G78KgWABQH4TpZyLaOO+ucw6c7YUf88dv8mHyXzSvrTXlZOQ2AvGM5xLMDaJRpfg0beZQDSmIq02hWsQ4ZG/XNfIO5Lu/qgRYKvT9GvWz138yQkXSLO5Xz5d2rDjdq04tinCs/7JPXLzK7mvvxrtDVXcFoNbRXdJZSq9c2S//X7U/TrtWRenQZgCc/THID/d71tDGQC/uhaC+gRg6mLBvDF+FUKAAu6wtxDouD5ACRrEzU+nealWZ3T67Qc2tfQlnyL3UquXfo/t5m5EhG1ZZSrqZ7MEFYtTNskXnv0k/rnaZVsW4NlG8p70v5qaC24Wn1yOF/+nVqtdc+uwQ7vVcel6NK1FTdvAnBsO4l1pvFpANaBwy0A25cQfauFns/if/tgyTOTuTXYxe0W0NRFA/huBGAAWBCA76Fhb1qG9Gqtbwul49M6VVuatU5Ukih1J3AoknEtr9oAJE2a1r/K+yJHhapDdJ13vuZSp7G+O8/vtfbOF+qfp/lj7eHcmktXv21/NbTdwhpj5biz5V9t8lw0rs+Beq4UladF9HJVinWZ1TQH4Jq17bN+gaWFf+dsj3SuRfO/3ZE8oMuj99Xct+s/MLcGe7mz1EWScNSt+8RgAF+FAAwAC5oM3cMaXI31Vfl0PVbXkqctzUpc1OVfyX1t32xKblRH+xR7C6z2F4e+lCuf5lx8m6Bb5g3Ase/UtfR7of7ZzH2w9F6L37a/GlojaM3wlwJAax+thdzbABzkcl0BtocYC7GaxpN28bKk37pc6Wq2JHHdI12TLYC38cTaEVo/nb/eS/XXtxCAsZf/2ushL4n8WRDANyEAA8CCAPwj2/G7HoA0NgNPq/3Afd3Ua8rVpVGvk5BqbL2Rp6SFwfPwXqPHt/Bqvabk4CUi9n3AS+vpbf2zaclWg6urvl6pf576rCTdIpy3i1p90m8tlufbo1ffulTXEJ3Xu90EYH3IrM27dOk4lTYxSe+/lmxNwvTrcq3Su5Tqg37aFtNi/O8+WBJaWIjDLnihA3BkBGAAWPC+8EcX65/X12rX47Y0W1uXKcm6MVXNsF7zqnaD0ibN2jx5xNe+lNtWevsW4jZtaBnz2+htR6qcTuufp94Hy2YpaROs6UpWlAArJ6XLsScH2PAkO8AWsW01OrQKbd0vnKqVcI/F4WnMRmpnqHuP29qy5uEU+yAop+2veob3Xj9NeZqnQP1XnqVDG/ZCNzUAR8YrIAAseF/4o3X9s7XCWo8ask97XXPQ7bQuFK0qDvPEXd1m7Xw9qQDupc4tJfZWz9ryWBdfN8u82lnqYv3z1Ptg6UQlib8Suq8F4HbPsWSX4/pvHcs6dq2jijvq/mU5sFgG18vnYUtGq7Jr385rG5s1Rbe1aKsM16/J+f5VOKffmZCmeejxf/W0YoM6drFjNzU3++Nh51e5MzucLgA0vKAAwIK3WbdZhbNlVCtIXm+1tXLoaa5n1nZWbTFW823Ive9TnrthreJrXxcNYdyDkPAcJci6vLl/Da95GpuHT8g9uNgj9/kMJLvbrJtyS83asXm+AyvM7p/UagOENQD7qL2r5qhuAXharT2P4cCTZdrQlpHbWrRNitKNya61m279o4MrWhbey73/rw8W5QnYxV6bydevlneG24uHXcy3vBQDeBxeXwBgwbuu2ywo9gXSNvV3vdDaWzTb0mwLezG3tlVW/9z2vFpFsy631my3Wuqfc7Z7sAuTz7YIPGKq7T2WkLmuiz6Rkq7YphRSuRaAJVTrTOKSXOltsCTJr5t4CQvAbZtucjrMaL5crol6t+Np0gun29q1ZtrQ7ko+8t6Wx/WhXFv+bYOP5HsiEb2kYp2t/ysA8+TELnbZTH5nZP3xMPuUAAzgmXh9AYCudzzCdTrapy0C29TfXq48s9JojaxRd/Am19KvZGCXR8cnLVFOOgBpBNRe/9yWWceOYi0ZTpoRR/nx1PKnbj+uab7NVgqlTeXNPtRNpu2yzk+qoe3Rrb4P6G3jhddHya3lEm367HTj7bh8BOC5LdccgEcDal3adjYK2AJwK7aOPbG3SUi6+TgUq6O+uIx99ZtPJMAedtlMvksAHh9fTMUA8CC8xABARwC+bdQ/2yKwtaRaH6DZr1X2Fh+1z3PUg3NMvQ1yo0XQKUqMlKssdmqYTH1D77hDfQPctsnWXCw7W8V1r22+9P5Yl11Drt7p4qq/HIBzyPo4KVkAlo8tyW8O00ico1ZT6wPmk8uj3m0LvO1BbfW4BXJ73F4WLtG5thXgNv1IC6fbl5Os+1ZbVV4NgfrZjvs2cXC71NLvtQJ88UJ2/wJ4KF5WAKAjY9xmC7zWBEsi3Zh4ZOxTbeLssoRQrfuV72dxVcNoz4vW+akG3QncZvtGKw22HlHjDpeiaFs1jZof7XjdZFvctv1Vo2vLWbNxLVX3Hp//gsvtxErvF20bdDdNvJavNGr/audOooIG4ND7Y1lw1Y5fk18mMPl5KnLQKmzdBpx0vVdv0IqeJVRbDXZKVhN9byNohgBjL7tEykcH4B9vAgC/xmsKAHQE4Nus/rnV7eoqaJ9XNLNP5ToJwDbmVz6t3vX6ZzumpVj71KLj+pJxh73xs/We8r7tqNXy45hLiH0Y7/np6THVabrWn6M239okW5strGwiUduge7FZtI56irp1eROA9QyjntVoTd0D8LwnWRd124lITs65BeAcg5/bRre22PKp7W229HvnM44AjF3Iq9wudS4PDcC/OAAA7scLCgB0ZIwbbHnW9rVa6fJm7VRLo3OVuGubbC3bJZd0B+yq/jnlXHL/VO8tFD14rn+2O+z1z3arttiqgbb9dUJC7MU1WxtBZDN45aHt+JOdvbo7uCzlzE7XaSXeXvhSS0na+iq19eqT35I9AJey3oOsS81tSXnqIVy/HrmHkloALpp4+58A2ghiTR/ey8n4tnJ8ZxhJzV2HAtc9ogX0tUv+ctj/HgAA9+MFBQA6AvANtmXX6p+nOQ+vD9C9r3HqIbNt+pXE6dpq55j3q5Nvqx/rxqnNIrKNsNv657E7ti2uapqOWtgs14XiNj2rphaTtauz3LMF77QNwDWmk7fQbRyxa82otkKIOaQadWDw6dtuLXROTrtVp6V0WQN5LRZkx2lE/cr0e+VsX/n85Ug6tpLv6IvtjL7zjf0ujYuAvf6SQgAG8Ll4QQGAjgB8jUVf2wBsqdJ25I4D+vqwrxr0cq3ahVk7OGv981zjq8k21PX6bS6tntn3gmqrRtZl3liWXNjKlX0NvZGWc7YVeXOGy6bfVjCtLaZyXJc3R5dPfrZzAN4WU+csJ68bd3VKsB6xvrIH4JYgRojQRd/aT7ivPJeUspyyrgeH2hNvz/Op5XC9++R93+Z8Tx8sPd/7G0YDV+z1l5Rn7gG+disA+B1eUACgk1hDAL4o9t27cTRMHj2cjfaqahOGtBbYawesWLyTd9puqRXWvJdOYrMVDNtA4DFRyZZ5T6Yc6ZKs0/bRXgcgtex5soSltcqpDxAe24ZdWh3WYud6sTcnL/cpZ3gSKlufKo3xSTK4z/oVxfWtdExx1gA8CrSnMUWpTUbShJ/ly9ZeWfongLZgbom3J+ZStE21Rn+9vDX5uqsP1i6de4Edn0gXBxr97rA7pyUBwC54TQGAjm2W11g6Hcu/07TtsTzqn3UZvW367Tts59xpU380Vc6Lt3aJzQiyteVp1D+fptVcYo+L+vilj95dn57XlldLIJe3+Ln0LcF2Dz6t/7KhObZ42+h78gMPwRZm9SFqD8DrhKwL3UW7V817ftuNLAAHreVObc5R0uLu1N626/dC68DTvBLetiKXpIG53WJaLybfIHdGAMbf7Rsmr00qujbW6L/OijFIAB6ElxUA6AjAF0nqs6rjJU9eGoCkoTXpiF8rBu71z3Pu7OOR6lL/bJcULRD2o6Da2lxtphxJFrUJut7pbfvo3fHorf3VSSRuVcU+zjl5TD+aaRl2yRaAl4dqodYKkjUAlzbBtw3sHexyW2E+CcCpjSZue4P7ynOv3NZ4r4XWab5BrRqA5Sv2XsO0n9Z3dQMxAH93vq0dAA6I10EA6OgzdFGfbzQti7ebXbhygNYMS9ytreVVm+ir1cUubOqf1+vGWjLsQ83JtQ28S/1zCOv6577x2Gu3KOuw1ScPjUdv439PdgW3RdWQek6W2647LfewXaucsS9xCZ+6Z7n0lDrp5mS5Ti5bPyNsZVhXmOvSvVm/P6mv5Oqy9lh5btuMewDOyw2Sz7rs7ZwuU/dJTz/8CPYaXYODo80BAEwEYAAYCMAXSf7UFLdqXnVe/6xBLib5pzdDnoIbs39skda1bbdjDdkqgq1xVButNM1dsjZtkaO1r2rDhOdpu8uGYb25K9vBSG1RNeU+oXc9/cgKue3gEnVGUc+VbSfu6DzdA7CETp/XNQE9e2s18pxaa41Zvz/2JWkADvPKc2tYravI9rSa67pTKBaAxyP+uCZHbsEuaHMAABMBGAAGAvA5K2/WDb3z74vz+mctJJaA2tov11ymS/XPyWdJg2OdVrtczaXOo7eWpkEfN/XP2kg5alctWzK12Nm33VqQjieNtVSLp1qNrC2ZT6YfjV7TQuf0tk7NduYpLY9sOdk7J5H2QgBu/a763UrSlseQ+2ypW7totQCsQVoL6l3Kvu80nts951T7ZOBWcd3GNv3QCJoAjF3wEgcAEwEYAAYmzZyzdGrNnMYlm/pnXdr0WrcsMVUSoNU/ay30/M3UwBhO6p/1JvMBLVB6WyW2bDkOs+3HfZZSS6i2jDpqsLWPVL0wFljuq8gpyGn55Q2/DXMah9QYNNi7PqR30ybaAnBKdR08ewBugbV/1X0iUt/0q3N9daOzt5ZdEoBjiRKRNQDPra37pt/WAsv6YP0YgNmdjl3QSxwAJgIwAAwE4HO2QXfT/3ldb+xzW/5t+VU7YFlmdmVd/9xG+7p1/XNv9TyHUgm0PhZdEz7tB6V31XK1BsSWMC0rWubUhduY1vuBFyHoDtsSXCjLZadbl6vu3NWttSUVXaM+DZgtG7vzDlUaodsuZb1h0NvYkrjGCl0x1gDcT6kU3WZcg+4K9mW0e+7l0630Oekw4p8nIbFwh13QSxwAJgIwAAwsj2xYA6qxR3c6q3/WUuE29VcrgV2yjKaDcF1c1z9Lsl3XP69bPVtLKt1jnKMtI69PQLcfh7RquBxSq0nuK8NREq1fZ9pFS7QSgGMryZ7Oln+nFoD1q9N12gu7cNv83ssBuKaofyxx/a8l2pW6lVLbmrbOZLIHqrXoVmjN9voXAVvqtQDc2lnrl+PzmAZ8Q9Apy5e+TOButIAGAMNLIQB0BOCNUfw89s1u6p91ddei3TytUzNzKBoF5++klgqf1T9Pfql/lgCZa7HmUutHt7Ddl39NStbS2UYBa3q89lushVqN4lc6V6uikTj60gYYnX3tmmn9eYtm7b+VWktq30crabvp2vJrW57WReD5rFoptE5SksSr35P5vvR7Uqp+OS7YCvvtAMwzE3/HTnIAMARgAOiIGRu9OLllVLtk3X55ajHWol1yfaCurui6pZLZqp3X9c+6k9eXdf2zfBBi7SurK3pX8XRttGi01rXWSffZyrXb9lenR/aRvGejg6f55LTVc6oXw2eQJHspAGuDq6CB1aqmdfE2LW2fdfySBODav0s16l8QUtVl3vFV2/dNByPZSbYAfHvIEc9M/B0BGAAMARgAOuoD12z0kf3bLrHC43GAta2yAUUSg205VLfyujhWVK3aeV3/vG713EfythXRUH2peX0Cbfyv37aHatOAdf9s6zZ1of3VfHKSLaNGar3/C8u/dldF//9iAI65ddyaTpela9XzzLF3mXb6BSxzidsuX+2DVXqTavn6ew1567BlQ5Lt21JSsbHJVhl9+6nHMxN/x05yADD8TgWAjpixpvOHpjAy6rTKq/0Aq1l2LvtkGVLrn1fjfydLj2G1hThfqn/OmqJ1YtBcaN3vSiJoqNvTapuNtXly8es0viGpW7JlKckGJl1Y/m2NpyQA6wbdS3eTtLmXxtaTvtQh9KDe2nHZctoSgG3+cFoF4Bh9dX1Ocu2dn6cxBardoaVo+fe1/mtVz+HmAjFwB3kW0eQPACYCMAAMBODB2l/Jv6/VP9vuX92g24YA28JSzGU9/ndb/yzvvkPqXaVrHfXPtiZcclwHWl01vbheVfRBvRYZu8vtr9oDa2G2DxKAvRY4h+1CcRv8q48i2b3mi+kyxphaKl1mFLXmz7ojumqO1xZc8/n2hD/PE172Hidr8+x65I+6pbnPTfLZ7nN04bqWTeS7EG5vEQbuQAtoADC82wMAxTrbmsQ8nc07LV2mNl2UNVZa+6ucx98NJMYu439LkXRZWgfkGrztc00uaatn3Tos1+nycil9GbnWk3Vad+3PEVY27KqVZ188eS26Dho+q2TzemmhuIVaiZ0+h1jTxQeKugFZu2Ppsbm2uUrJCsI1fXsJssHLZS3o9tFQbcewrnLnuW1YzhaAnQ3Y0rHCyWqe9S8FVftgybdRvlFy0LVJSGzdxC74Ax8AGF4NAUARgNds+fda/bOlVsu6kjZHOuvjf3X8T8ouStzV5lgup6g7hPv4XyspbiOCJKBq3yjfF43HUm2vf74W+tqqqStXOmBJ5rQBRaXUtlScS90cYEu3WntcNKmeFDmvHkR7V7c9vcnHqNuRdVNxKG1fsc0ZbiN/rU9YXxtvk5NG8y05h5h7AE6t65XcMMc2+1f3L+s/+v3x+cYkpNRc+V4Ad6GOAAAGAjAAKK3P5Q1iMxZ+1zXP61poLWx2vQJ4tL+SXClxt3gdAtyWeDXjxRB1gm22Rd+28bXdqs3OndeE50VjKw22x9JZSvXK6XldNV16L59faxfrnzQ0fG5Li1tEjy3Qy8PJf0O4UH7svQ9Z1651WG/KcmStbd27tJOc13Ll+zO6henNLBjH3D+tNaXWBMsyrJVe1yrnr6OA5SI5oP2ZIF3pxTXRuwh7oI4AAAYCMAAoAvCgk32mtOn5vK5/1vDmXMk1pzpWzX3oo25rW3HVoCfpd+7/bIvGXhd99VpbT9bSYr+M/7Xpvhq/S7j2Xl0eMfgsiVEe7uJkI60ozu3jWlvRp47hXX1tumm5zo267Gtse3u392QBuAQnSbsFVaub7tOJ2zSnOEpK2/JtWFaG4/ytqzXHHoBDa/tse481PMv3sH2REuVtN/W1+gMCMP6OZxEADARgAFCUCA6WdSXOjS5TJ7XQcdJOTk5XOKPXtCmxLoaqDbFOx/9K2Bux2VZcqw8SgHU1uIY6aX8qa6Nlx1g7ZTnwcvurabK1U+2x7IIOQTr7FVZ8HDFSTl7jq55rma+ee1/Nj2mR27pSbb8JztkatX0tY5XYbmLFzKN2WjtjTdpxOic/6ZrvsjpdYquRzprY7X6tAtz6Y7cDss1SvrZDkyHA+Dt5cSMAA4AhAAOAokTQjKx7sf5ZZ9mm5G2pVxd0W7oLufd2Ph3/u/R/tgXh1OufY9YNtTb7p40APm3+fL1ZTwmaLbXjVNDV13VVtlov/9r5By3I9mm+qKXY3nbL7rBNSMpL4+ouaWesHoDt+JGQewBua7kjFVvxtnxzLACvV6dL0h3CocwBeL6NfkO8fuG1fTdsY/PFRtAEYPwdzyIAGAjAAKAIwEY3zU5lXf9sW4L7ZtfqdYFTG0tJ2NOpvykECYoaS1clvDbsd9Q/L+2vWsLTRdeqbZ905TPEcbldpTOE8qUza/2i+oShUvpS8CoAF7/UJI9CZTlbn1t4lYjakva64NnWnNtq7vI42har+GDji52zGUXDMvXXuU0AtpOU1O3j0tG6Rg3A8vXK90r/oDDfxla/tSDcVoNbO+iLva5o3ou/41kEAAMviACgCMCTZb+W3Db1z7YwqxuDvbZ3tu+TNrVqrbA0bWZdE7YQqHG3TQMeJcq2IDxSpl2uTZhd6yDd9tPaVbqqXMLllSpJpG7eitx6Ry3ThqYLy7+ajXXWkZPDdLOxtacu02arrWb7UteF06F6fSq086/eRV/WN1kCcAhWAT6NvD3NnbfCMjKqZm2RVbTGWUcfT3O9tQVguYdJb1FLKvIVndfg69Qnogv+hv0dALDGr1UAUAybmVpuPBls28jHthpsMVLSWi5VZwm5vnKaU2v+NFWbY6TFvc7nsmwb7vtd28qnZUULydoly6bstohp7a8utLbS67It/y6jmFxYpg2dLv9ab2f9qA1b0oViHy12rhably+5lB5xe2/qUuyvIXqV9xJNNwm0h9vwj703MJDbVpZ2GQpCQQj9Xd4AACAASURBVCgIBZkgFISCUHi6qhogyJld7bFlW9J2P/1+EpdDgpw5e+djVVfnml2zvQRzWaPLuG4gHeOYonxgt0X7nCnMMVGzwbil+pqDFdlsUX+/4uleVFRU1F4BwFFRUVGoAOBldb70zCn/+h45t6NoQJGx39JSDTKFlWqpRT6W6JH+ZDDrccm/2i5IVkeuOmCBuIy/2tn7quOQ/Lt+hDbglTWFkb1tvXv7bgJgjBq+Kc1XrfW04fR7zg9DGrkZj7em4KprLTp+5fiiLUxrLbUaAA/japqujV/7oVPU0c7ptxYAjwJRvWTkaWtk8IPPA12i/n7ZR+ib/3KLioqK2isAOCoqKgoVY0KUwHyyEVc6MIj09D5eDTLCaN/eDNjOyWq1tVG2FKtejePaqAsIPe1pyr8eOrUhMQg5AQshNReXkW8rs6/vRzU2XliOQncv3ci9G0sv7fSSf1m2NiR1ZQ1FOl/N1Thdh/c4NfxlnrCk1AyA7SeQcO+6sceD2cdlusFF/vOI4G1McmpMBWMvNAiZM4TXJKQloRsDaw7w6GNP25pXWQOAo/5mxS+3qKioqL0CgKOioqJQ3/w7ouRfkecag6TWX+0wcvHhQxlDa1vp2o6O2dktjFclQ7ljWakd8yaeqrUYOHpcTIkRuwk9srkU23YFSqkI3gf7cW9g3HtOwGkQ5hYf9RCQRwXDU4h9M+z39EQuiMnoUl5XUUrOsmpTws03AHZE77gPovirMZiKNnRyermxkjHwRKDD7J01M4p3A0qw7oztcLScPA3r4XcO7S7q71dEQEdFRUXtFQAcFRUVhfrmAExC5IwiKrI+1XZ6oUVrCF5OQ0OA9XXabtg+jBf37zjG6IulfdLvvLGKg1JE1n52+7kdE526DzvxOQXV0m4bz5mDZfidr+7ffV6x9jk5zcj+30eNtOB8gOyFmSRuNAYDgHsRAO8fDQdgXGcW2mM61Ol6r+2PcU2tLEW9qtk3NR/IROgFCB9jzSUGDKPd95nR9c0/llE/pSJHLSoqKmqv+J0YFRUVhUKr5tsxrN+gJP9Ks5V4Kz3z8kLnYX+U24zs4smYho07cBr4GbHp5T5eaJN/hdN47UurK16Sm/a81NRT5A3IhBv54YvWWOCeU/EoZqfQvWYL7vHxewsFGPnL/bqoTGoeGPtUB6Y0lXJzJq+I7N16jVMLsu1FxuS8M2owri2VgcQs24iroKDsSjBhGOBMGNZN3rW6hEHEHyw9KuoLxYy3d8FyUVFRUd+1AoCjoqKiUN8ZgJfVWWnGDy+0elMxDejwFGi/T63ttmTJvxjqQ+qz4xia7vKvcFoRWa8LSLWtnVeAFuzEpe9jdfcCZ7bisuqG617dDcafADCCqE4Eb61RwBk9v66YOQBzFPAOwMukPSQ+w8ONYchuvTa4TVkdy7qfbDGG/ZnLqMuNjbuVHIY1ORma9v3hQJhXo/5mRY5aVFRU1KMCgKOioqJQ35Y0JP+u+KtFv8tOXPKACfmAmbkdZXmJEU21/R8Ro8El/9o/62iGzKulFXOAzmNB8r4A/ciYcsnF7jHGJsi/R7+Trcp+au9Z9VHAtxgqFYVWg0mEK5f22kiroUeg047JTkhprr5eAbAPNyKE7wRxdSlzNDH6mdkr7NBObVcJ1RKla3M9/MDYpOPUmbTAfT5Trhq5tC81zKtRf7OijTwqKirqUfF/WaOioqJQ3xaABb2IYBKtzQAql2E7kCwnNu4yBdr1yVr3qbyjS8S87MGlDgPgPpp20HbjvTEF230BUnGNtCECQyQlYR5Hs4O8lX+p7vY2MIuXvcHP7OjmKdP2n9I4MfiOAN6obOBr4GkALEhVQ+8Ysox6N7IBcu17C7G3RnMN0MPbuC2AAGzLBgAPPkdoSAWDzJtGHunU5F8BcK47DEsKXufqtshx3CLBoqL+z4o28qioqKhHBQBHRUVFob4nAEt9daV3G/nr4ieh1OhX5ueW6uolNB7e85YRoUyQWxvrUVtNOojQGrD50oto+OjKMKlWBmlBqeYPvZd/OZjI/lNyUdPyMyKLOrMouI2aermFYLWm0URaG9zL5fIeGy0o3Nm7kW09ue4vN2zG8nD6sp4I3AzYx2F4fpCjDXprZYg1AVhzm7C8Al8qwrR5aLQZ2x3DQwZe3OgagHQkvAjTlQODo/5SwSdR6zf85RYVFRX1UQUAR0VFRaG+p9dU8i+0WZLwAlrHOcq/4LGDjayHx00p13iXf3XzVmp0q0Nx0NJFFSstdfexgDIgz17/zJBPjYoNShEN9Vb+ZZ9uY5J0L8UDt3b5l022a/Cv8Woa+ULvgtjnxZMCYOU5LwBWz6QDMBqBJ/nb4WrtJUHIPTmKad4TXSN/jtvlI5r0io4O4Vxwl6Cxs69Y6da1jJGgD8P+jGnAxsq45wdUa+rWeAqQ+OQh+VylqKj/pw40sdfA4KioqKhV3/ELX1RUVNRrfUMAFvSK09bMnrXd/qLWXwwoysc+bUjzhxZzIiOaL10guuKsVqpWqQN230ehI/bY7Zn2d7wWcjySq47+Ev4846oEt6NeDcxem8FYZK32YF98RxfuQy4G+R9pDfuV9Lp+5KFWx6l0LPDwyqnmKCadBYOjevX23pxL6qtB2s5pl1IaybZnD9YCVVe7rxjSVJPdX2PpclRmUJ+tD4VIE6qHneOwK6iHbwkMjvparQhofYrgzS8lMDgqKuqb17f7whcVFRX1tr4hAC/5d/1F29WvK1ewUdlIGGN7yb+lKOLYj0KV2L5RLx8yrM7TMWyoZhsh8z6jr/BCNMf2Wz4tApUPOITRGNzL09g8Vd0Ft6NfZmYdU6fu1xJcyBXf9pIekjL8yCPxGv2Ye2gQe3fh3MZpZ/fyPqhJAAw0rjWPeS3IB8MzApEqlPDOK8oD05TtJ8eBe9Jz4sWuscCYFsxJSDr/NVDq7ANXneGI7uGIjvpqPSKgdwz+tqH3UVFRUd/uC19UVFTUa63co+9TknmV+ay/rB9J9sz04o4jGf3C0yuepOa5y79rgu5KgUKy8fbdGiJzy2/msGQmTr+IUZ62Va886v0lCna+WpEHMpyrorbm3KPTe4Sv19mh7HJqKo9jCqgxZDhfo4D30CDAZ8234U/3gcOYRWxHyDjE9VAAii1mFJdWQK3woOIOQ9NuJeGZAC6kDUyWMvb2scAHMsAOaN83evdZTXat57CFIXirHNCs2fz8clujoq6Si/6xMTA4Kirqm1cAcFRUVNR3BGCpvoKrh/xrf8R1mHtUU88b/pWSSr2E2dYyGXbJv62OR2AyxvCO9ORcfC0/Lsn0WpaBZFWA1rP7txQNEF5wK28zEpp7VfayfjB3vAoG72KgXB9ELVkY6A4+7nMU8DU1OA8CcEpX/DWLrbwNQjRNy3303NAYbPfGcNf4tNl/7KCD45DschkyLaEYjw+OOf7XbkzB8e2NUK/1ofcjjb2x2dD3rSNaDvZwREd9VJ9EQC8MtgoMjoqK+lYVABwVFRV1dcp9k/K22BO+XDHw+hEaYivCjdH6y+wr+4vzZO/iOqflMdSvi4lEJfVWZIfev0vTe1yfrbzGfkfaM7fWdqNPRDfTd30D4H6lJS+4FSIi5HlqqueUgR+FH9e8k/l5NQvT8m3HmElXngcONzbkZVxsVgBYn5ZmnPpInSlV+WDcsx0fFwsGRogXULYVzXMi+gJWc8VLAM92uwoU3dVajB2Mu6m0p4P26H4Jv1rwckTfMLhFY3DUh/XDcPsx3FZgFdOSoqKivkkFAEdFRUUBgJ949EcX+lEnAK/JvadinEc2ygWCHogmbjld7buGgZnZTphvW+COPuCu7PgKDRIeMDu31ewL+bfh5T5Qd5XhYku35KrzIteVO301984f3czPp08hhoW7e+aWLM2vahY4v+YdjI061xsOTG1ppHwBMCcV2UY9IICOdtidwG6FivbRQby1CYBPTVGWGGtAjsHCR1J3sFZogA1lt/sN7DVjplTph52qdj1EsAvGnS8IG0M/cPUO6r3jdzmib43B0xEdjcFRj/pitIEwOOPxU/uHVxQVFRX131cAcFRUVNT3AmAJiXLa7u2spzTVfHgMFbASwOYthOxr9ZBiA7Z8TUVS9y+G5SZOAVIU8gAQ6uU6i9MsDMv5ePxfn0muGIdbGZWMnlyqrxvUrqBmKdj+UqP1kXwY75tBS3PWcU2LBfaIrJMAzHFLZPUxDuDk0TGnl0yPdl7MJMJoX4DvkNK74FxdxLIor8sRAF+TolrdAbh2O0D2e8lnB3577U8ryh5bK9yFXz8aHdHShzUx2OB9OaKjMThKNfhh/q9XERUVFfXLVfxmjIqKinpmpf7ZJZoSkWp6rbZD/m1oXpUCXGpDTLGYriGYSnHQEnh7KiI07/6FEjobCTkMFzvnjLm+LD8RbceGi09I29KtmBpdKYFyWu+E2l2zXYFb2FSAgopifn2IIfYuA43F6hG+ZGSYrTOM3Cd4FElXaRAljwI9NkPebRBgdVsUziwN1sC4DpqTB+zKrY/VUz3vcrIbabvZqbEbnc95FGRcnaTgnGmyNpjGvcXc44RnCicfMthf2Hf9Xvhd79cbR3Q0BkfN+la/1qKioqK+XgHAUVFRUd/omyIod4qHD/k3M5xJIUyA1A4s642URkH46v6tVbFYQLhROFu37vBpUAlEPKYavDRbo8rRnsONONRIQm/tbpbGNrbgytu8m58h2EpZJf2ecnSX/qaJe4BLO0zGHERkpFog2LZmi++1FEZTYRrvAdMydij8F1qCjYtPN2+XjJ7eBJIdYlE0FffDzlkhQDv5Sxl2+DTOzwBUOKJHagbShriDo4BPUnACykovh8zcKfoe0Kax8WAWdD8ewu9yRHs38gbGFwbnyxEdGPyda59oHRUVFRW1KgA4Kioq6hsBsPp+Zdbd5V/7CwcfZRmcke6ce01Nbmfon9Azs1RWTypWWDGEYdBaWcOPOmjT9gf5jssUjXlFve8jlFAzikpC71I4DXuNRQGKGtV7XOZnh/ZJv6cSrZXcvOf9kH5rw2CkNhBGXSiPtgMJ18m4l2UoPDITp4/DrjclfhIYUm243Dpk8JIR9wULNC/a7htDt9z2XI6WmoddtbPpL3hkkKvuMza2hAWcVIAx+Ohgd6aBb1PD85qKhGbgXpWD1XJZfLtQVsd/BeO1GxqDS15ScDQGf9v6JAL6r9Ux62/uFsbsqKio/7bid1BUVFTU+2mZf16JnWR+fsi/pRHJSHi9doRZHY6mI7ksDIwtw3jPY7EaXtN7pYLrDmTuk45HSzUd1NVeOu7ZVxNiRcEr1wonzQVoTdV0H+oLRXc0icb4t3qDEVldmUMlmzIwUEBf2dBrQClFGpm4RrwpX8OU6NamMdlA1nb3RyE4RKnwg5NXbRHQw0drFdBq16KlQoPNndHPhyzQCsRq8FK3vSsYcu7gKOAOWMVjA7u0XlIedsfsXCcDwBBHjpNh3jIeQPAUr8S7wHgPMFu7RWNw1PmFCOj/q3Zq/SLcvt3tKwgdFRUV9Y9W/A6KioqK+hYALOL1WGMFIy/5dxQ3OR9CWWRQiX57hhpp2FZGVR6VO34L4pCx7955C+pq4Lw0rjAqb+3NLZcLudX/Knc0I6n2XCuBMaTmDm/zTr9aqtiVxMuWYgNgcGurpYh4W6b4WWuXUTsRYo/Z0ixdWkxvJDv4LMB2TgUDiQjArQ+NLMKQ5HplRCvn2S4RVm3SRcOVDUi5w2kTiG6HOZpeq3tqB4N0bDsRgPGfOQqYlJ1h/85Qg+GqTgBgOKJ1pXer83pDJfzafz9zRCfekWgM/n6FOLefBMCvyPoR3H6+Rf8MAI6KivpvK34HRUVFRf18r+AvWLLsqsN2p00GOHvIE1DWyC2XFdEESMzDRx+RfgGtdrvKwfbajX4p4wLAqKGitXVQ+yU4FjsehVY4nHXo++Aiz7XawBiysS2muyaMzlt5kEsZg028FYpuSlg/NWdssf9/zzxRScsRjQsZqeVt2LN81a5ZQxhHuzO0bqPRw+h7JVpBrW2Qk9NGE7iBA3OSgM2Vic22HvqcB6KvGehs/x+7f+3v0I0r9i/JRwHbDeytkM2ZQn5w0RwFrPts74idEW3GH8xAemwxLF+7LUd0O9uAkl12R3Rg8Depn8iZPwWA198DgKOiov7bit9BUVFRUX8+AK/RRwKnPUwY1Eo8Eok12oMdD4GXQFPA1YHBRCC91kZOnOXTL0GVmi0adOdEIui+x5CJFyjbs1CwjMt+vNNv5The8e1cF+KSkc9cq4BzIJ/5wDGYdgyhF82zozWcF0OYkOMMnVkarPce02cNwBbEJv6ba1A+8ylj8miYXXRAuMX9KR7f5cKpYq5T0UbJ6a0XNU6rd7dBCc5EcBwQzxQ6qF6MinSqkWqBsJyGRh7hRPCTd7y8UNtWcphGIqFjmbd651u3Og+PMTs/0IfXOy59eOC6UzQGf5/6ubkGP0sB/mRjVFRU1L9W8TsoKioq6g8HYNEa4prYtrq32nrA1Wz9lStYKcTsmD3bQKIVkLgSdzk6SJOK8JI2XEo1Kt6H/RKJkcVE2MT84JOoqm7b03dhy2tBkHJ3N7IhpBiY4jMUVYVGgTbzIUszZFgajJdDWB5mwG1OLjKzbHkIdj5oSD6pWh8usU6F2pcCwqychGQELM1Zzm1EUvGAkJl7La0zKhptu0y6AvZjXjC6fhUZbS9F8nNnotjMzbabo+nHuSKvSwAM6O0clcx7aBiM41HQxg2ASxrIrR5goazW/JEjWru91YeXIxoYPB3R0Rj8B1cAcFRUVNRHFb+DoqKiosB0fzAAyxCr/ttHq+2ae2QQZnxYK6f9DrazDh8XLDs09skQS42agGaiX5mkGUDlHubel/R60hatlmAJrUDijtwsH3rUBkhPxma2YdvPE3tfET9lfJjZHUtB9Zz51K4359FSFaMCU+U8boWDe5MvSWZqnlqKq8GwselJPXQt8qS/Wn25BsC9Nr8/bBXGggvkYhygGEQWRWDDGi0v8SC4Nria99Qr6MmwcE+5u5FFjcd5MbXj1tnR4MJumbBb8XghwSueR/EcrKNKqdYKMVHpndVZGLx2W2OuXncDBjepy9EY/CfXz32oFwAcFRX1J1X8DoqKiooCAF9RS39WaW7Q0mYdU098QXb6TQZVpKYBYKs5t1IRwlSKOAoY3DOs0bV2ZjNjCHB2u/IoNOvaDp3SssCK37zZYIuzl1RL6hxHhAOVmoTTPNpBBCYNHvb3cqDldsBszHG2OHIBarYGy7X06lZmt3BKAmNkLGPhcPvake2FmvqrmzAa8B6tuUd1VWzLwUI8NW6DUSiM1R5/PfVhScdA6INXAbTuolkPuAKx2iYPx7p6j3NWu69ahfUSpV4BZSf2G5Pbsu1V4OeOO4ML58RgWwkjx971AM+Ny+q8UPZzR7T04WgM/uMLaeI/73daAHBUVNSfVPE7KCoqKuonf1n8pepyz1LndPNzBSJyNJD3o+KfTGDCUB9DUHTUHi0DI3NqoEdOMDL6RVuvxNiEv5Ots6YoCVntjxRG414wJ/kWOArSrEdBNBSyjmVsLkifMnSsqWCg7kH0JdnazsBadg4vRdcOb8d3zNa8ooF5uZlLVyJXpyQLjpXGy45fI2DZkrPkaNagwZuTfhkPBgUcrc4uLBNua8+yKz9eqyHDCHNOuFn6+EikFfHisUp2PZbi9KGpTmiuhtCLBw6G0Ea/uKV2l0qSsqyxwIgJYw6WnWIXflcP8EdWZ7+0zRH9kT5sn/ndER2NwX9S/XMzkD7a8nd2i4qKivo3K34HRUVFRf3kL4u/Ti3v60kpWBgMS+2R8mxJtW+jmDubu/RG+yM3MeDNGNOQLReALhk31UJ/cnX01VTbnuXXNeIVTgGhmd5EdRdULJbTARHgVJOPAuJPdTTYfUm2IEC6mo2YU/F+Wp1CaHr18VKPxYlG9vG5CRQ9Wr/kaJm0CYQaL6zX9sErarWzsxc8aZxudybPBmC+FgL3IGmv11IfhjJ8tt4ZEA1Ftp8zXaw1JGm3BAC2U6IFWoxNNzWanQdAVOFbRvjgah4WduWM+2BXkznd2B9G6CnGeBIvniCM+tbq/Nhtp+XPHdHRGPxn1M+FzADgqKioP6nid1BUVFTUnwnAC37YmIvgJTiTDQXR5wpqramMI4l+NQSY7mJSaWn6p+uD0IhdMVYmldGl74ks4yFuTUeT2Htwsi3w7jgff/BT9hhnzOxFu69OYX93rCWQ27JgOE6l5Mb8Z+AoRF0/d1NmlWTkboDZKRcfA5FXfqaMYUm167XaE7bk4+zoj8Qmo1SPlRq4Pw3J0li/BhermbmVA328JEUA8DJON+fMPpB5xWwsPg7IXdOFa4aszUlIePRg12Mntf/iCmiHtn1HbwzT7j5ymJnVrSAHC+SvHCwjX9tNGu8Ajf+/xHuLwuJrd0d0NAb/eWX/+8jXgLKfU28HGv213QKAo6Ki/tuK30FRUVFR+EL2hwHw3vqr4KtOcmvlEJrKhGx0eQNU8m3OFYJvQzIyKc14NWvCkPMVbNGV7bvgV/Cv4rNS670B8wYZOw07zijIdrI/HohVmhK2Sh5uZGae85KOM3uPbwSOIb6uFWPZE2hzdZ82A6tBsMDFXJcZG6cc3hVsW+Q6BoFz/WozXl3EUFMNGo6kmUbex0uG9pDqnH0UMN3YGGvUlFbtqrCxOjTeOvuH7ZWpqAdYqdFM80K2tjqAIYAz1Br9x2MKzrayMijRd9fVmWJ9bsS7O6LtHfko/PnYorD0Gficlr0xOBzRv3/93AjoVces1+1f2e2j/aOioqL+5YrfQVFRUVF/2hcyTSQSeqmzVGZj2pHZ70zNV2HOMtkaywmQtANaU5snLCscS+5fdAOLW6GognglKVMcTcBrtgprI8TlWuRbNpAEblFGLon6ZzPOY2h0G5KgMYG3AueM/S4wTlg9fwxnsei3M/XqgGXZKRHhWEaGqWGMUM1QZ6UYU40dvWGME2YuJR+MJHJd/Apjsv3U4Jufhemv1s2EeE5ORidzHUq5xgijTjN2zpi3lLFnGWiHHtn7hxEzZjeEG6EkN+yNl+ehpw0aR2xnA6Zm7GZ3HtJdSrq0a8jwvQf4I0f0Hv78OfGu+UwrJFyvbb3BRX4GBv/G9WfPdYuKior6m/VHfeeLioqK+mv1JwCwcZV95eXIHHfGtiyNlCKuwqmAo+TZY4YyFc09AnpN7JF+ipFCmBqUfGjtMRtiKxpuacfNzNDKmZZm90HXqh5Xo1A7MHRc21KnEGqIWKh/jqRoK629j2Gr9RNNKd6uBnLxga0GroDJNKD0Zu9G1loxjYmZXgrcEqi7SjwVY7h7NXyoJc0W2vXhAR5VjjRcwKDc3DxqK7l42xHx3BZv2z/RUz2SHNHAVM4/UswYrjQdSqoCvWONtbCNGkhZgNAAzp6UgwXQ7UiNRp4X5gV3wfbBKC/gMd8+EPdb4uUN+yHxro0IM+NG+5x80REdjcG/V/3Zc92ioqKi/mb9/t/5oqKiov5eYWLr7w7AnDqLP8lBtJL4NGsHQi/6XYsIFplMxE1DHZFPyWX5im3nRu3V0Gs15aI3lZFUw17GP8A2IZNhcGmZ05D87GwwXiytvGKoypSCjYwPnlrRVnYo8RiWMbCwkfKV0ozpuQNRUjA4p9YHzg15FDbrVLoUXDQFV/zIWFWadoXVeDIcnwEktjsbgnok9eHSNF7M2cWQlKFB00fdeUs54QnbzoNtumeFwlx3MMYqWxNCQ6dtGQI4jdMUyx2AdS3wRI+D4VjIJ7O7aTukgrcKfuw5vlgMj95sJkKrDRhNnR1vX+t3jZcG9TraeSdevF/vMPgt8b7udmHw7ONetBz169cfGWoQFRUV9bPqN//OFxUVFfW3CxKg9LvfqEpRFtTJ8TtoSq21FMdOHx1ECbRQwu20M4/W5S6+cLQURDCl6Xw+nO7syBBLWxJGLg+zDkUZmVIlYQl/E7TWKS+z29gAtXaXgnEOTkvSwF4cejiZI2O6XoFStnSdcW9z1eQkiKinRh/ROJ0GifVqYFZItUPvAadxITDrbkC7brgnSyWGwj0xGNOA+WJEc9mR2aIMJk1NWdfYxwA4ew4WLup0qzPeCHqkNbLoVNMzd8N7oalIHbOLXZkHFlPRzQfJ+0RqdNkc0bng5uSuFWKRpSyN126Xbs5FrcLgDWXTuBHvD3uA18bliPbG4BPzkcMR/XvVb/9ELyoqKuqfrPgVGRUV9d3rNwDgFzHHEMhbbqdTt5WJgUV+5JRTQ5csRWAInhpcZHwHHIV72QOTee1GbOr1NfaTBuuABECFSVoblSZtDCww86ZW7gZBM5+jD1dHZZBWIy000gnGdooC9zXHDiGFC8bsnIWOdig7RR2Yj6sty69bZ2bVAmO7K6UWNhonewlG7BpjA7VdNNYVIZCZIi1iuXoSM1/dv/MPloeIrObDh46yIrIcgI85Rni+Fp3Mc1qSUrIw5WjMYxYMNM41KxJbG9WGjcsf87BcIWC+ZR0Hkm9jxzUHL602YF2z3NR4+cAkp3PPeaZN+m8Srz4MiK1Wp3M0Bv9u1TFE+idHQEdFRUX9SRUAHBUV9d3rN/i+WIoP5jGQG/0UzHRIvuhLnXQFKEpN5l5vea1Ge0b3TcZgDZ3Vy+GkZXMsEHqlPZ2ed4XYZCYhe8hwnUzIWbWGumBgtOPCeNzlDtZuCW5kl4JHEezh8ExOFhirfVdjgYwYDbE8auuBpFPKzvyRoSQnBhf1zS5+MzbTRrQin/aKgkFEJ63F98PVDPzLpShbC2BchqK2KBfbXTL2rXId4wnC4ZI4lGflYE/1+GTELk+XpvncnxpcAMy3DKcbyR3RdJhXZo8xJPuQ6tuyB445h+euEC84qOmC1n+lD893H482kK3Vi56PqCdZwvJTH1ZiVn/2AO/E64FnMwrrc0e0gGaINAAAIABJREFU8Dgag3/N+ocioKOioqL+mAoAjoqK+u71ywGwRgYRdXxDJTkmtI/aX+AiPrfcqUm/nvbE1CjNBypsgoX9dgzkQafuHbYZNGVEdDmfibtL9TWy9Y2G3ARjNuiOZYc2xJIXGg3A/RAYIxB5FFd9NUS3zhDpRHF4gTF7lD1dWag2EAFtF4gG41G9MZixWMpPPmQ8hou7KsV6tbna1fnspZODndI4OcCosbc5VyQC4aLL8CHGeCww8M/Du4AHRXFqujePtNRgqLg99VZqBYtmTPhFnpam8moisQB1MO4L/Jhm/pYafm09vUm3xxLSaJw8rHfPGNnubSN7D5wGCrYeCKw2YKE4NhYnXui0JN5eLri9xgW3a+PSeLHxTrzS2N8S7/5aEa/GCOthh7ZEY/AvWEa/tYZEHxUVFfVhBQBHRUV99/qPBRPNAtqqtCpMEqOKWIzx1AtK9bJn9ugK19RwK423JLSULqMvBF7DTp2HovGao7s2yqvsZFubraW0jgPm/rkd2kXdUpU+hZf35nZo0PFlhwZELTv0koI5Cgge4Dl0F2Q1xybppMZpS65UG7AxJNBrHIt4He0Yv2z7aNJvna3Rdi9BodWFTSicxqx2JGZo2ZVqplHBivKRwKIY3guFNivyGlZqiahUXJGwRQCW9q6kMfmZEdOVQcWj8eZo8HLK3ip8CsrhrBZUOyrfkRvt0Q1g2mllVwIWU8W6tGgfnjQPqHZfhXjjvJNaxcDYSMsA7qTaoRvDzL5MvOd0RNsN1MZoDP7FK2YgRUVFRX1eAcBRUVHfvf5bABYR7ZIvCDGRK+GYJWudyrDKSlQWdSTmNRN+MgbzGBzmpGZW8ZFBIKjJ/1H9GLVcuNsH22iHUZGbdg3Mjs7eVTa6umCJRCj7Ro1k5bsdGqJxBl5CND4uWr6Itwyd2u3QuXufMFVcXF7zHuZlw6aG3JcxWxsRxXUingr6MHlMSIzoaczNdUcup+ziuEjWmv5eUGcCuBeirG3sTALLac4Znn7rI8/07IOPAFrJQGl4ltVKfSojmuHYChgzjmUYF5RbXB+fSoBReXNwbgTyzonECHnmiCO7sdSKe6VuzF5f8bYmKsGRgN50nEMPPpyQKY8v7/RyRNulGd/aRhxwPU2gTaA1eNQveZ+XsXvIdTPBt8O902+J93JTs+vbG4OHu6YDg3+digjoqKioqM8rADgqKuq7178KwJy4czt7n7oqJV/jKO5lTMTWVmqD6JrNWXHK+iP5l1HKtOmmrkk/SHXu89AGW2WOMqKZ2TefED8X8YqObRkf2aGFwTVxJC9NybsdGleUMF4INM4ZuRB4GRl9NQBznz0Zy9Yy4HeGEReTkI6GdeYJb7ADl7USrRmImpnaVS+7r6vBbEa+8AzqZPfZufL3TilY7vHasYaOVWakZLesoLCjHJzv1PMa/kOSZaszZiRJV9cUJfsDczUai51RGVvVUnd3usKf7Y1UVpbGSmmZdgfx0OH0Acu4M+mUr5t4jJZgOdU7tPS+moF1w2VgbiRcOaJ1NxSyDem4XbdI9+fVEV1B7f55053EFhLv0of1tGUR707LV2MwTpEYSl6iMfhXqIiAjoqKivq84rdkVFTUd6/K+nfO1cC7m4WVBTycDNwSXc2EMQQ6N+i6JbFjFulG1fOigJGabeQUhaG+MPTOgyupmOOJ1mAhOGZdCabw+yq0bnZoATNjn7PSoSEyJ48sNgC87NClDCq0BsY4+Xy5M/BqAKbzeQfj2obagCGTLuLlebWYfSwwMKx3Oo7zAuNrz5byhFbIsIbJBfqwNHNjM3XMLjqF9tuLfgRozXiQAMYdSGD2AcI0PefW5NYW89suImSkZxn89qy/i+FByI0dznP4cIOLObeOqUgQjTnfCJzMY2oSkmzbvbD3m2+77Otg40LtN2e1Lq/8rQW3uyN6bTSwx8Vt2WZ29xrHTY3RF9xK9VUf767xQkweH+rDouWdeKMx+Ncp+5/Ir55pHxUVFfVfVwBwVFTUd69/DoChZE4jrkpSpOzNEEURf4TwpZQQg8QJRoY0iAz2kbnj+iMhzggNgm2DOpjr5I3R3ed8zDTlLdYLFtlySr2Ee7n5xrdCKyiIJudFnpel+RhG4RR43VItcVh2aDCS2o+3PmGnLw0KnnbobuclVB9sApYNuMrRfYw9lXpfod22QpswA6jrWt65ZUGDYw3JBlF4WpuV4ey+cT4+aBlUmhCMVfWm4ENwcEQU3y/FMldkXrFb2K50KuHqwU4cvgT9d1qjQdbpkCwMY/aBjuXVRZyU80wLM3XkgtSrmhTOjGPnOZGYx4HsbDcWlAwYPrkSD6B+OMbxfmA88eWI7m5gVtSZ9GEYp88kKdju5O6I1m3Z9WEp5ys4+hPi1ZaF0N4/PK7+4ah/uSICOioqKuqHFQAcFRX13eunhaYCF27HARLUKhjoGBaEMuJF1hKJbzBiScNvQZXAlCLmAf9klwRv2GNUyawpdXh6FJFP9fFWz9I4JQfnBCqTi6fku3GjVpjqRK+5XTOB3Q6tRt+VDr17pIfnCEsaleUY7bKTlldk9ErGwsFav6RgAiJQjXOefEhSqQ8mB6lqLLBhWGbA8hi7rO23fwwp5Cf16oYhuW15pFsv4j1ooThTh0LcL1SzQ6fq2VoYMMUHCK1C90YsFhuYr05pQ86qELKkK8WNVXC0vb89+RhhtQ1rFLCQ2I5dk0b74o6JPA/2E9tfS1qdxnhWQOI9qL1LlPb2Y1ionXgvFm354YjObGQGjU5aVqJYsZvXrlFJIlhEg3ExO9za/VGM1kPjXcHRD0f0qd7jgQ9FGdEY/B9UREBHRUVF/bACgKOior57/azQ1N4g7c6m09EweseoJstO6wm/lVQqKTjRE3s05SGJFmC4rQUtozQbu105D9pXgYboth0YStQwHKm1AnW1Us2U4LdgGZN+mGi1vg8jXqtWsNmmEX2YO7XboReLyg59OEjPqGNMGAIWNjQDI27qZLvv6hN+TcZi0JYnY9nl5KFFLzDemVxSsBK/NCvIwQ8aatnRXTiHG3VcFODwRsO4gA0IV4tkaoVU40cYqgyNFNljCTZShze2CivtGUttW3twamJgPlywa0bAFcYpjaS8bvEqBg9rpvGcKryecUigFiTbWyCXMlzQuFNNF4tl5aZpwDoChG7mfO0q/TI/7/fNWBdXVK+75FFhA651nKvd5FyYn+uVeoWNw0dAfUK80oe1/z4qKRqD//2KCOioqKioH1YAcFRU1Hevv/iVsTWkDXHGzMnIXsi7UHYh7BrBZAQ4G6E2aL+NEVPj2DUx0iM69uSz5e6uMmqsDtVLTcBxbAMEktwExtcfipBG4Hv46+jDpxNRwYOzF3pwK3Vy49ZG+3U7NFTE+0aM7UmnHNuabTtKdsUYaV3uat5TtZYderSu5GoERBVPq0Z28gAQ7slYtn6gVyIbT4DHNN1a9t00ixiia7rljUH1zYXjoniH2SkJsVcjkQcGHB30mbfe7IfQzNkkLAc73NcV4cm4obSRK87MrdHMbebUKKzQbwJypdFDDNG/H56VzTZqqNkHjNZyuWu8Mzo37copCJfcNPdohZwh4nu1AesWGVi2/IBbfMAMbrVxg1vYnvNTNFZjMHaWceAsC4MlBS99GE3U9Qm374mX/mcflaQnCz0ag/+9igjoqKioqB9WAHBUVNR3r68AMJtvN+cwQKWzAxftl2ziTTIf2g9Gq0ob3r/xgyUYloQO1DV2lS3BjVOPoPcybFlZwQ0uW8IvTMN9nheU1aapdVsPpF2ZXcW6VXNxT8dpEfXqKxYBwhjMi1JsLAjznoCFrtjpfDaWQ5csLvXaeLNDZ46j7VmSJlBoNgAvO7Tm1j6SsaSsopnWeFmR0fJI70w+c7AAt+nKgj63BmDtRoaGNC0d9cb5vNuJduhzULLu7iTHM4TkQ5Wh4Q/ceLtDIl4pmcaQCuJWgJYbpyFkV6VGu2OZ7cdYDx8EuJ48dHV8rtEcgO34msykFC57IZ4psL0b84QPKM/QezkYiXOKu0TmdXvt7UMDcH/jiDbs3B3RnVfaWlYu2rm1+9qRV2KWEBq4O+pb4rXTPTbehgPP6cHLJo0Dsjc6b27zqH+uIgI6Kioq6ocVvyijoqK+ewHBrtlBKMX87Fv66ALXPHEXucXVIIGTaeuddYeStY7HN3413MJCjKAk76XEWNo222X7kKOYAUgzyVkYY38oO4L2SudaGHpkJ+sY6NMgKJYlIANkm3KOE2OSDsdavaQ1YzMDbIE6pD/bARzbDVQgMyaEIq9lg+Fmay5DmaYdurQ1algtwew6Ta+qL7zGmx06nWnZocFmg67mnODmFhgvj3TPy54tuK21Cwt3uBUtQy9tEBvxlh25tiKBd39eAD2TkVHIte6nIA3O3mOKw5q0nKAt45EB5XY1D58u/GPalGBP/l60WzMCGifNlK+BtgOdxqU2WqAFrs7A/dBoX3jm2XVsb5BtvTqHS+HNKcxFYyh0xuMGfFQ1TFkfDI4aNqxFLFetu1CfNLzpYX4mLXuC2Xx3wLR0RAuDH3Iu/M+bI9rbqueoJO8f3oj3I0c01PVoDP6HKyKgo6Kior5SAcBRUVHfvV4BWFqWfVlXZpWI19AXiED+uL0e43aT5LJ9M/KTKPo9bc/TG5zRFZqFoLD/ClYNTYdUwJlZXKrDHplw/oC6Wjku1zTPgr5j0dHdBnlNSLqZpJEoDVNxPa6BOrLzijCZBGWX7124lDoXG6+NF+6O02GOec4rjPpSfd/ZoUvFS6BtQk4/rj5h2aEN7EeRFKwcrA6pGrC/D0nSegw28XBhJMnIWpJtvA0ZxrErsPRoOMssu8zV4GqYV+laBvINn6+78BhvFft11+OD+dbi6UVNjBArXffSe3179ujvw2ck+fTmyuZkDkPKnLSsLt+Txvje0OXsmvOBKCy7iJYbVOV8nBvbi0tf4fbVEb1E41dHNFLDMucnjbHkXFx+d7jFTKzpiEaS1l0KVtC0/fcj0dgDtHo0Bv9TFRHQUVFRUV+pAOCoqKjvXvDyslZKM4gXumMxCnywsaKGdjs0ivAj1r1Jx1RFIcqQhpaTmfNvM9tK0dWJjKNycawMt+jXNbJucCpjCDB6R1vh9NyV5ISDNOYe7Ws5qR7fZU/fXtOyN6+NtXmGc9lYvQjV7hZiIOumx14bNzt0q/RbKxoq1SVUSi209ed02aHXkGEcNjfh99UnXOXtviKjgYgniBnG8Jw+6lsmg2dNbPILH3V5pO2k0nJr7y5f892XdIZ72xl5ZWfL3habKP/j1IPCbye3n0Vm8vX0wV7eMxm48cODvDIHYE1CUhNvO4p6j5klRvs3ArIPDRmWMoxcaPrhNWM5c7IvXmgXmmmublWZYS2XRbyCW+PYV0d07dUfbajdl526Oohs0oJbdP82OKLxI33at3bfV7jFUKUNgyXRy1q/iFc3XH+RI5qpXtEY/POLKQQBwFFRUVE/qADgqKio7172lVEarwB4Ee8yiD7UVLHck3VPfP00jn06PI03ypxzw15TYTDwoGd1q2rgDabVJuYh0b2s7sr9D7KF62F0onlCC644IujZEmzAs+bf7CsHmbxsx1QhBDjlysE8hmK6LxPSbiqr4qnepEMXH6eEhtiuTCgIlQZ1ayaw4LbUmx26cjoR5iHlIVlSIViSgnubUrA2coqvTNrn7Kx+ZEETo4vvszUAo036SGrd1ha2EsM53DIKFw3irMDNIzvJb3/EnNC9j5PvXad/mzq8Ap8003iguRczlua4Y+A0lwQ1mZ3S6vJly+YhRRcC7+Hubs++ouqrJxFg4IJ7a6/HRs7HwnnxqeguN9MRjQ9PL50Xq8lYSD4b8JPv+rCk3Y9E49eNS+Pd5dzliJZNern6NZ1Yveg+v3p2F2sj2q2jMfhnV0RAR0VFRX2lAoCjoqK+exmiKrH57U8NLlrJrySpHsuHFDz6ZXvek6tO0oiBysJgtRlrbu2gKicMrrL1HlvL7jxIna3CknAJO/jrGhWcZrNwZUwyorXsnLmrd1YHw2TgAVA/70qR+BBScEZkEZgE9wWq80PyRQgWnc9vJhVVZzYQJcU9tKEeSXgvPVZ0tNuhz+WRPk4R4LnmIal5uJ7OwBxZJH81IHA7tRqAV1ZWY7iYW8f3uCxmQadzOocThipDOE2Z6uqQhqn3XaN68ZItTgxGcVK6ZghLB0agVuEoXU5R0p6Y87yuotOwrvcIfFlyY6ozwqrhf8YReKMYgebQq7lHLhozBwspanl0ZYQfh5Ynaf2cPcBXv+47uEVi1mwMXqKx3YQhJblfWVYIDy+Y8qyN3uU7uu1v5CyuXnArm7Tg9jI/c+TyY+OO0GgMbtEY/NMqIqCjoqKivlIBwFFRUVE+xEga4Nsf9+wNsde2KV2CYe+NteplVevj2mxANtjqm5AFVLx9lJbm1Q/pGKzBv+oFXVxLfjVOgX76IuHCU50xvKfnajuKeAXDaiu1A5bZ9Xn5cmmJte2ibFwRrcZo1yWT6Fpqr7vTWDGzbofeIPbc7NAYkNuzVEE7AdTOPTLajcqXkmw3+NCE3AP9vb20Ky5rSsFIxuJYJfd+H7csaHD1ZPI6yjIkt5QfcVkQfk+keuG1OeviObEZY3Ux8vf0Nxo4Wi5WFBYuuFUOs9i4CZz1Ts0wLWngmjMERzQS0Ibf/wyMBPM2Ii0gNCvsCjB80CB/+Nwj5GjBpAC+0TTgwvnATQ9KCqC0Nc7vRR7WteAy8qsj2rXc3pYjWuZnQOnmiF5w60Ou0E19cWxnHPbDEa2g6d0RjccQLz3Ai4eXLIzxyy0ag39CRQR0VFRU1FcqfldGRUVFeTkFvZVQFpu+tAQbtT2JVBN+U1KK0oXB5ElGa7l1GQnAJ/yv6RGUVT0BCy2Wol/xKztGV1DwAmMNIoLFNx0PNkYKVHm2BEs9zgXSoiGWDNWCwaMkhUtpApD4B/bjCofzuM8f+sQObX+gE1ZXkiF43qcH15mMpUZcTGkaiIeCQfjIVwgWuesWGa1e3DQeXcouBSdahdXlm71V+OGRZsJZ1qQiKeEYYdUqArbzJH+ko50NTOkGY3uJdHu0Z1c840DH7HnFQZeGNx0Sbm+rMdgxuOUtqMzeqS7Fm4Lywc7qUxboernbXRKXU1p/kSysLbdHJBX+gNrRWI70tcn88jnDD49R1DeO3X3Okr4VB/02I/pVScaN0cE/gFtvKqYsLLiVur4gWX3IgmQOTU4rGj0ag/9CjflwKioqKirq84rflVFRUVFXKQzpIzu0sdDeEOsbMTpoOmM3mVG2Zyb8pssRPTyh6mr9ZR6xErCUfeVNwisLmtNkL+4mvPYKzHNsnUNxEImVqqdJbxiv+be7VOvb1RJMZl6nw/Da4u2yaDk+y2Jg+HDvdujlfL55s7WxwRvsQquxUgaUaibwskOPNQ/pwHBjDijCTnYKn4dUXX6se2Q0e3ERkoxU5ProUsapQce8Ia2rVRjv0eqJZaOvuq/LNAOv5klbko7vGnpd45oTDeLTN66hTUwj0xunnm2DOmIdBv+K+vSecg7VYXy9pFEIwrmrX1qZYQbkoNw2AfhA+zQjuaaMz6FKiXZpPCOgVg8SznWO0Woe8oyBWuXhiLaf7hyr5XUy8gW3HaLx4Mv1ho4562h3REu5FeW+ou8nPIxZyncFGIbq2TaMxuASjcF/sSICOioqKuqLFQAcFRUVdasxxicMDIKVy3QnTMOJnnyM0PZC8MiR0MDJ1l9JW2JjF+iaG54BduwSVXLwAie4iIerwfdjw/Zc0ibtToQVFD1WeNKD/crGpyYkbdsPztoBatKGTbQiwXK00mvw8ikt967u+kY7jHYWn/e2PNIPOzSYDuNz2QebIMgacF7zkJYJeRRFRrfeehaQZgxJyu0hBYsnEXl1uMkZtN9dCl4eaQVhkTyvuG/QFyTesgZB6eVIje7zMYfeFUVJzwlJrvFmBlwN4O5lGx4cdzSwYmNOWamVv2X4Kvs3/8Y7mZsGAiv7yv6u6GyXyo2ZIbJrIHDVQxBq19R4S8F9bvAnQ6SdnxAQ5vCNyxG9fM6A2+SicZ0xzrghdERfojE/CciIrreNdkx8/sebHmAXe6e9+SuOaIRj1WgM/r8rIqCjoqKivlgBwFFRUVHPUkvwR3Zo+DNrfbCuvuXnwvGvmxSsQT4GKujxnI5o27nXPlsr58HKNf8X422YvSRyUIqSHdwwSXNUvWtXtudXA3bb8p/uoO5EdI+KHduEJLkoXQ0usAQD8IzRhEMaaFTbg4FvdugNGn2g0VHXTOCrrffVDn0eEH7ReYsHAwuMn3Zo48xW1d7sXP0SypU4Wjkxc3l/OwTwUpt9kcDNeqSjz8Vg1DCfO9Sa93eZKz90o6hEux/eoa7O1l8Ku5CF27H3D9eCpQIOe2KuGsYp6R4q4lu6LswAjMsS4srtzNFS3IHYjynErT5c0AameMZBwvfPD83PyxEtffgrjujd5/zJRvvQPDhW07CXaPwKt5KCdccemVg7OcsRXZh5lrY8uahPKiKgo6Kior5YAcBRUVFR7+szBn5rKh4+hPbRcSrRWK2/0O7qdHiOsXZjp67/CzyX0wqCptrISTYnEoNdbxx0PCNRGDNVIeE+MrrI3a4Gq1V4bv98QlKdbYQi85lonHJGKBNZ+FAQNOzQ9wlJNxBVG+1JojNu3+KyXPXlSOGnHXqyH06s0KY1D8nQKWeJmWoeVmLWHpe1HNoOsYZnua8u2bVOQ155pBlvhtZUuyZDQ60c7cGnvykaSrTG52JkER4W+HQiAbDsu2gV5syjFQR9i8saHYlcdVLfOBCbPJBlpfxkEGNLPnuJYVx2oYJeXBGPl+koL9D4GbTGTmsxsAYpjWJv+Sjz8YQtFDlYJza2zREtOsXt3czPaYY8743BrxsXx+KGtNoLPqLwXU+x1zi82oX0W7vvmpn0UQ/w655IF7MrLymxFz0ag39YEQEdFRUV9cUKAI6Kior6sD6fkLQaeneVFbOK6PK9ScFT58W4mYwZOEozhhc3lUsy9Tm4yF6yHeyPzziaQCIbM9zJRGJYZY0Y2UQqptpDdEXp3iqc0g7bztgv1+UTku5xWfY6Ozq4rGRfovTAAZEOBLvt/7BDA1Y7M7EGsrV2b7bboTf3MnYWzB3nq0ca85AKnhNgUNMJY7ZPHuZ6XB+eDm13blPKhr84l0frMtZTwISd99+jrcGJ915fA+jWXB6ngH/1cucbPWoj0pHJqxp2JauzYtIaG2iXzomfdhBvbzKgA0rtTREA44ZXxHqv7CsBsLqCGRbNZR64gatnGBhcwL16LIIlwB0NWwJgdTqiNZxp+ZxlfpZB3UG0D9+4+Zzx/vAjjdDpcduID+R9T6H+2mg3Te2+skkvXXefn+Szl5S8RZzWa492UBefj42iPqg1zPy/XkhUVFTUr14BwFFRUVGf1ectwbQbPwVVVyNXr+z6kRtH+4q/qlBGk+dIT4pWYzCQjcQFxljTfgl7O3U7WRt2UWh1GJ562nlu04Nrdbzmevra/lCNxrRDb4FekIJLQwcse0edOacH+yM7NPKQK6YEpdIRHzayZvz4QTY7tIKstDM0xuLRx9dqz34Zpwtc0B1I1yWQ1jkoWMidOpquDbfkMNcyrnXy1JoVbBsNKcehzmNfPEYhzV5fNgvXpfra9mXH7VRoB7ulj+3/mEIKhm/b3wj/M2CNXki/0rM4dCn7FChmQisWC2HUdEQDeusgSXsbsGThkhCOBQE5TYmbjuiRJgZ3Pi45JvECUhOl4OmIxvMCsHGH4/9yRK8hwOc2WFgTlXYpeOzjguf/BHbzs23EZKZ+ycvLEY3nLFNJlvVAGPxDR3Q0Bn9SBsB6YBcYHBUVFfV5BQBHRUVF/aB+0BK8hNY7BGZNlU3XTBrfWz7e7n2/0Ouq67dGC+uYAg0gH6fvzDMx5gr64M3We3uJ/RDqXtUoWkzu2VuFl9maBuOa6+M4WvzbyUmiekCR4pGZuOPEmNrDZryLsRh1m7EROCRb7/RO73ZoESYlRFwFtO0jryZedVnLDo0hwyWd0wXtkdH9istyL7oszQnE7kvKTaf2O8ONRp35vni8lg8L+qbPG9JlZHZPIV2noxQM1bZt4WcZdxtwWB3qOAyJoc/lWiHeoKnwA+nB1JhvnCqeNSjBix2/xo4N3vijLEc00JdDhP0WUQHWDp4OnQif2xMKjTjW7b0c0X3O+8Wn3DVtXOnceJmfp89537iQdfM4lH0jvNAN/cmvcKsosrXRxwiPd3uedTmiozH4be0R0AuDoyU4Kioq6m0FAEdFRUV9qcTAH32nvPKldggEhuWbBntOV/EWf9XbhcFqgPQjSHqlG/amel3tuceTXYW3pSorWJpwYgqxIQ8oXd+StQYkBg/OdnpznDUhab/GxJG+bifeGNhbgjeb8bns0Kl7wlPn1mzQt8Vl7XZoW3amcjtoFG9u8d2VW+1p6Fc6IByDcw/YtnGZ8khXl6MTByxhetBBOGdpHtLukcYxE2hh7/WlsHka7PrDi+2K0OLrYumpHCxpoer11eX00uTT7ulKUbZlywyOxuCatUKNboY+PACodjQEQUOAhsRtB0FX8OGfAYxAOlwKxqSowwVwb5nmH83U0vYpBbf17oySvUfXPnT96Yi2M753RM+MaF2IvNO0Rtxs0vjcygSe025+to2aDrUGAqs9Wz3Agu19o6RmvhlZDfBS0ZcjGm9XNAbf63UGkm0RBvdrhFpUVFRUFCoAOCoqKuqr9XlLMLH0CYHqluSQI1pYhyPB4rQrCYupQpoJvHuYFwY/T7sw+CX+CsjNIUANJteZxYSRO1XC5qlBPkvc+2hC0st2BVB5VyrpbvEh+kZfJiS5yZkNt75MgBPn+OS+7+x7cnowkpkxTgi0bH97REZ7XBac6Wgh1XCgPS6H+DP0AAAgAElEQVRrUL5GmvLZsFEi/F2dxnqGv1PFqDsD6lZ6mZpRteeVp73yrsa8jTMHCzOrRl1dwYr4ZvpU6zNKmnp7XmHRiAsfB3pwmYMF/3NjXzBzsDwLiuSpjl+5oNUALN1b6KscLA/l0ripWpDyxUcM7oju1X0K0xGNLGiXgi9HtLfj3h3Rzpl97I5odfM+bNLL54w7wFZkzf7VRluDERn++4iD1v3cZia51Iyg7bJLwRizPB3RiGJLRziiV30UAd3wSUiBwVFRUVF7BQBHRUVF/R+FttBSPmsJnphxzm+cmvvaOmKvcp1Ecd4YeKKoY61Cji+imJJteg4w2oK47kvqYPXZrjmwltKbY/DKm17HVZtrurGrH35NSJqXAzOqYWr1mbregDp3viYkrdzpaXIGtqUikXoYdDIy+o0d2oCdg3BxTkqp9d30YFp+ixTR3fm8pgejVXqUaXKu3ns8r04eaam+xoOVl4CbVimubr2+yJFK2YFw5l15kNhRpIFTDfdWYQitJGHR2tp44mMDIVQNwIn5WBlThPNySmPBJN41VdhuGyg3YVW4qEN2+aZmYGViyVWuaUnbVCSmqN0d0bpAbYFAOByDe3Fnu/uQ5V7m5+Ra2+Zz/spGNQb3hkcDa/avksYf7b7a4WGT9pytd47oqzG4RGOw12fmlNnBUWBGGP/qsqKioqJ+yQoAjoqKivr/6vOWYCukRN3H88olCx8ytdwyQ3ThOJ27XVLweUraXRgsfylAg/lKG41qPaemwp73H3jr72ZjHitRWVLtYmABTB/lJdBLrxK3XAIs9Vvjo1yaQ9c2+/fthCQZkqXWum483gxDOgXYtEN7FjRVVnc+ryCr1qQHKjAZ6D7zpT0uC9bz7qldA5FRdteg8m6rggMZ0mop6J49Vi4X0Ktl6LGackQIX/qwJiHpCK16rnWHzjt1+IHOarTvtrRvtF2QsC3SpjR6y8ES0gOJPQcLKvEJIoVxYFq7L0c0Lc1LDZYUvIzQHryF7Ku2HNF4+lHpx9YFHYh0tkXCn4wZwX3NEMa97dDAYc6fjuiV0nxzRPfZjjs8wdz2XD5nPi6iI7rffM4VI6Hwkb66ee2mMZFrbZQC33qRTTqxj1obdfzdER2NwT+cgRQYHBUVFbUqADgqKirqr5S+TX4iBZd88w+LGzl4BqQAf+fsgVzgu0vBGEVDtoRYV5O+31eGERGLx067Y80lPq5O1LX9NaBr9cde9mkxcOtv5huz9glJRilQkjNUXle8dd51c2YC1jo+LKyjyMd7rXPrK36ou1fm85GutK0zL9GYWdi047YmEbgvxVsTiTleyK4IXccHSLUa4m0mbeVgEREJwPfuZY7zzfKHi1r1QCHtWWUdumvhvCIpw+u2IwEtHxcA69FA76tVeMndnMaEHmBpm8ruSrVJAT7pkQbEMm+5ZaJ647vPHCw9JlD81QLgNR/Ygbb2JQUDR5v6z90jbW+r2np3R7SR59KHlyMaUdIzRmuZnzHKiB9s0CnZWhs9MQsfl6aTXjFap0KiXcvdN+om2J+10W3YI9dN6dUThPkQBI7ohG7x/G0bg4/jS1/nfvjwLioqKuo7VABwVFRU1F+sH7cE3+OvlofTwAm8eQxMmB2HfbmnnHpNA5ZKart14UDy+at1RgpLDvXdhp8RLJq2ucGzrmVsDLx8xW5L1omzI5/bmx/fkrcJSRK0tUtldNMjketmh57JUgbMDKm+TSRyXp3DkFZEMzYiyGmsI+/LRvLz6QZydPlmQ8COEKx57e6mLtjYCjASDZFnW6sSQp9Uz+pRloasE3F+b73aa+d2XHjJy7uOn0qgPm8KvNF4L+jy3TdCU80uz2qRaeZguRp8Ah1rzsrBQhQWiJjdzkJKzFzCWzYMu3lRMoFLDd5FYPt7USzZgQHOpF08O7DtCsfCu2nLyV1bWnPd9VwjjqcjWgrtwxGtyKtlfnYb8/jMJg0M3jZKbH+bEY0W5Y/Nzw+b9NUYnL9pY7Dyrv7rVURFRUX9NhUAHBUVFfXXS4rKJ98+ofcanpFgBVH61o7v9GTAUl0JNNi6AGxJwWUmXc3G3aUcHlSHNSR2YbC38q79N/R6uLLPlaicYPn9sCX43liI7cXxJjHcS5cOdk5jNxifL3Zo8luD4be1y8y8iPHsaxhSaYcApnL+MHKutlFSY8ZlOXMykaoXyLBrUPCi9OtJgbKgWxtzHlLLblGGeEhWNKhekU4a7YsFYyrwdV2SstHNyz0hO/NHAOM5W+jklOCTsnCZ8qm/rRmW4NUVDJV1y8HyGc4Fr5xK+GFrBfR2538EZSfy5OFsvFzQj05gpKDdpWA9U+Do4+RSMPXdwh7jVm0xWRuXIxojrji59+Ts6OWIvrTWcd84sjYuR/TKc5YPXMHUy+d8clQSwuFoI18bF/Hujmje/CQpeDmi11Oh3RH9rTD4NQI6KioqKuqTCgCOioqK+lv1Q1eh4HCXgpFFNP2l0/AMRyhop27xUQ4RHmsMVqz+l4XBNFGPyzg9ZkvwarW9poOyZfdlshEQscBFDGY7J3xDcD5bb29agueEJBmzNRJ4vq6/TmaSHRo9qEQpDUAGkBtr3acZLzv0GoZkOx/sKhXgvY3LQjpUaoaSQN0VgjU1XpcHCcbeqyx9WB5pKrTHcdgtxwoNiieZe163zpWvDC0FQS9EXN3LsuOuVmFYr49DB4F1OfkbZ4TpmWGkQTtZ6WXlYDkDd8i5Yj91BfMC/CGFASsgFkOakFA9mMimgcBLBBYbY1xwPpSV1abEzeFHHo6FPTM/Tpz5tBzRc9nnzPR6cUSnfNNy7xu1J5baphS842hr4urbxjGQQd3TvnFwKhKk+JH3jRyqfKymgLVRg6+XI/r7NAZ/FAEdFRUVFfW2AoCjoqKifkJ1pCx91hKsDKNLDKQlWPDgVNvoMB5oyAQ5zNJPL31VBHn3VO/HEQYvjXEx8/lJS/BsuAVlaXt1EdX3f5G4MWCpctQQbLSu46oLdPcSr1MY6pVUla6s5uGho1OlfGRoSd2FTVfqbgO75pakDy9g9rgspDtBYpfAe8naZdgRpPGyJxfdrb20K3xbEmvqB0zCVUfWAgxigZ3Dnx10TrfS6CNhvN8EO0N1ZZ5C8GrkhmiMHKy+Gln9HcEDjtmuOdgVjHtyXDlYlfHUKwcLUDeQfYUzVzdU208zzcklO+4uF/TifCdhjUGCx7tK+MVGScF8xyWtq6uZ+jzuMHqS+yGb8e6Irn3TcnNdiVlLtn2Yn/eJX4+NmIdUL/MzMt5Ol4J3R/RJyV0NwCv76mS3sO1ZELd9bdSApasxuHyXxmAMC/sZAHzM+su7ffSj46X+/mqjoqKi/nLF76CoqKion1Owgn7MwOe0Qy/rr/DVvu6fU+wFu0JL44TYzb3pEHHMebY77E6WfsXgzmnAbr3OFze+bQnmOKB+Q1cqdVj0+HBCEobZsqX59sOcX/uQ3XjcGC/sU4f8so3iHuuBHfo8pBtD3S2nJidBi9xSrJSEJHqXDVggjcsZxaPChl8S5HOFRS1Fna7s3HOCZMiE6hXZNZDyBeBU/lantZghWG/yrii3LtjWApAmdtRHDpZhGyny2GdE4b1Oba3HLgqxUiDWWpqbog3ifDjwgGZrS50fFbywFZiiwdsYj3uypze7ZL1agssbRzRjyfqSghHvTGf+koLRoDvDsbRRsc8PR/QKx/L/GWCXpMgr7bnMDGoqXt5pNpyn1crrvD2nIq2Ny/ws3N3Nz2836oVYRLf7RlP0H+2I/imhVjuUfgKon+z2xR9FRUVF/ecVv5KioqKiflp9yQ5N1pS1eH1fZyfl1CY7pFENQNrdmy79LlYc48bAzR3ROpTsz62PZb5dyqcv492gpsZuZbywjmu5nCRcJjTeLoevAiZRH73uwgFx8jVNWvKmTLluze4dem0vu8P5FAMP6He2J2iZFz/mpF/tLC3bN9LF/YzLolZa5eguio+GLZgvxd3L3aAQHcBQZlPfL5DjlI45Gagp7GoweMkWJnVad8h+1BqSsQDGe4t15qWdzwcBPgxpbmTzcJOQu3K5bNXKwbI1qD8cg4760A20c2W9zwO4Cy96rxoIDPU3cegR3dGrJVhAq27qa2bSfVzwdff48bgyotnWiyulEm60K0e0bNKKg1YOVpqDjjwjuty0XD1S4Qjou8+5u1n9YX42/n/4nJWSvW/0buEx26fP7d2fndDLEd1G/SMx+O/j5VvN9v/a7fMjBABHRUX9UhW/kqKioqJ+cv2QgSXnthkmLCFLrHspu1TSQDube3O91r3KC5pbw19nV7A7Qikr0m88Dzp3nofqjynBYmA/RZ2n5GsvW/ULA9siEcJ8tOsnrbFJ+Ck1n+S9VIuw1gOwjE4xfak9LOIAPzp1oe4q3Zorryu/2liOnuqTbmTg2wrBFuev+cMJJI8oqQPngmhcmm0svTAuGXe4aALQhtBagDcwH3l1U+MSWl5LVQ6W9vcpR7oHdSje2R8qrEcPjQLpFNWpQleFVBciPY7GRl/1ta7hwKX6RYFajyEIF2YDAtkODXGYNmYZnn0S1czEkiVe290RvcKx+AThFo413xHNSVreY63czr4M3rsjGh+/cQm8fvyer3CsufHhnVYc9DVD+PTg69eNeEzQ3fyMMOr50dXk5NeNagxO+B/lH9gYjNFafzsC+h8F4KDfqKioX63it1JUVFTUz6/x6YSkc2mr5EPNUEWj4+wX9eyrTQp+OKJ9h9r831NRROQwoPaQTXcBsg/UlX162qFfW4LFDMYGVXLgjtnkutcW4nN1aVbIkms0sVhL/bSv9uZChCvAkO7S9j49eJKhfiyJ0j29ugOKy0JbcbEf6uDMly5OfeuGqCWY/a6YzHzUK0m7Vg+YmoQPj3TKy6O+9MxMzXldBa53y7tCXzEXxtFQWa3CeLrQTyF6A9+ntSTbBwQoQZUDje1a1sLs5uBK7V4e6BAeKwv6KLbaBZxoe0Zi16FmY8Rr8ZJ6Kj7yl4ZnAD9CtXw80ur4ZbM5aXk5olvFROXliC7ZFeY551kYrMQsjQt+OKJboyyfjke6FYiaGxdFYyOf/qipeMm2mookX/1jo6h1d0Tzg/Vikz7zvvHcGoN5G1MuGRjcb/v81vVTIqD/aQCO7t+oqKhfquI3UVRUVNQ/Uj9sCRZXorGVxld2rl6ge0nBVFNh1pVctr32ptNOh/OYLcEXRWiwcL7aj5Xz7CdaLcFShvnycTKBWQzc2uP4j5Zgl46hYDKMqm5tn+Qlpk9f9uY15tenEGsx6+xzTq9SpnUkMFpuj2FIMjkb3ungOEyh0lvzEpPHDBtTMpZz7NY7bVxZMj3SnDaEDlt2xuKfU14eTKtajc0Kgj4pO/uQJH6tVw6W28thFYaLGADMhCdBI+zvrawcLPS7FuKi3eeZJg2VNcMTnviAA5FO56E24PU0BOOEqk9FxoSnUvAuy+NtnMyAaGm8LdUlBQ/ivRzLtiCFY3lv+btwLIWx2VtoV+TjgqvjLjG46FAaTfSRwOuO6HzPwTovR7TafS8xmeHYD5s0+PnFEY2O4nfm508c0YBhTNsCCf8ZjcFGv5/8kvli/ZsKcDBwVFTUf17xaygqKirqn6oftgSfGnVrYMZwIGAh7dDz5ZcdGuD0Tgouu7w6A6LHbAlWqvA5AflImyf545ZgiWZrAdBp9Q17MbAmIW0txHKrulR8zKnCzaVRaYboO53HN6jznx/cnm7q7mTjlGbcFFRVA8eaHvFaux3as7Vao3Uc02B1sTAhU8jFRk3KLVXZYEB3I81cND1Ys4J9AXd5Wb2+NRV1xq6FDQ4AcmV1tBU6LX24+rzj+ZYx57sieOvYRzTjbdL7MuVl0B3bgDE4iZq5MTFc0A2iKfpguweJ8RKaUqwUhWU7gb6RY501HgnKMIVio027q0DZaQdA7jJDpCUaX47o7FIw4JQ3qna6npG5dSybsRzRa3jSwxF9+8TOj1wv6XWjKHppuVJ9e9ts0jzXkoIXRY85Kmk3Ng92C6+NyxH9tjH41qX8G9ZPmYH0T/cAf/1HUVFRUf9Cxe+gqKioqH+2xMD9Mgc/S1ppKV120BWLde7Kblew0BspuEldq9uxZOjl1305QtfRgA55TvcVoM6WYDDqZFq5SdcCpDGe6xBkwn3/5Z328y/b8xoTfJ4KRqb7t6cZ1pXVAisj7vY9XgiKcT/FM5yhjypeS82rm54MBqZ3Gl2+bSyIVVwWhNw+F5nFsW2J3jlDQ8ZMYHLVGrMkYEYQ1eTtxvu2hiTpmCrwHoOm0j0gGlOF7kHQAmbkSyude16C3QPxJ2LL9PRBgVPtUFgX8I/zgBHBxWZh6PP8UMlrDUN3zZK48XHiiCkBsNukGUxVKAirMdj+aW/YLgX7G4Ep0P0Kx1qaOfvGVzgWk8TG7oi2/4fF6+PRmFuWfbSvq753R/TKiMahpqF6Z+PGUcsPR7QmBu/UivlJ/CjpA/9qfl4bx3l1Vts6cB9/88bgnx4B/dGWz3cLAI6KivqNKn4HRUVFRf3j9ZWWYIBuHoVjSzXmdPt+P4ckNRhxX6Vg2ZidTDd117ancVfMXicqbROS9hZfSZF+/HpPn74m/152aA0l0ikgNasVVktfauc4NfLHqGYldWFnoHTbpd1zIqgIkLOheJiCO+StqjyvOMd7eg/q1Yfrq+58buBVB3uQPbRlDQpmq61hWx6anJTnYafZW2itBCwkQ1F1heLLZOa12pLhPS7V1nojWzxkgC573Bzjdtp+7FI2m6gh0goO0QZMO7GmAWvksr0XVzT0yRG4qeeZGo0orEzYYwJWL7hdfq/YPu2u8tLV+gvrcvErZeIVLQYFz1dWu/VyRCtfuhc8udDlG69jS76rvjVL9d3Jdgm8eo8ejmhcRZ9PZ/otI3ppuXJEl3FjYz3c0cOgtz7nTzaek43nnSx4CPLbNgb/FJgMAI6KivpWFb+DoqKiov6NQu7Pj7r1hKZ5zkq9ZqtOtzPYqlPhTMdOyNqjYqTPhLIJzep0fT3aPpF4Z+A+fDtSkVZAkdqV8xQDxz1JK/l2Nfee8jbn4Tbpzd48D/W0N3dqx5gPdBwPBpbDGc2rGgiM1mpPh17r1KE4uCg7288OZ10+Bhpx594473YeIcELCyXQXb7G4VRQ1xXBe3w2zQRacVPnnOfkQneHug+5tsMC7U8QeBUaPowo75r2LmhEIjMlW1K2moohIxuy9uZjluyNQEZXnf5thGMhDmvO9q2diVN6GlKB/QjbauxzPg6J4YPTlLV43QH08GIyLkb49lYUCu0zgGkxSGWTgqta069xwUseb80+pmWZn1107e6IhlF7XKqvXRo2tqf5+Vr/3RGNoVT7I56GNmFb7R6KvszP+oRdG2l+ftikdwx+dUTzMUD6HRuDu/6X9LfrX06BDgCOior6byt+B0VFRUX9S/WVlmC3Q1enVoh/G7VK2QU4lQseLsFKDPw6wWhrCda3/+toiE3KHos1s2RXS7CijJZ0nNXoqslJmw4MBkYC8CGhTwcRIIOBE8+8BdUOJkLt9mb/Dg/t+PCW4N0ePMqanOTnbQ0X1YvWv6dVIYc4t8cwpIO6ucYg4bs3V64MLeTSbnN95ZFGPPIcFCz+BDAntivzQiR3+5AkqqlS1PHjURWCJSO6PNhICOtXILadWm3AyruSdxrdvrMNWNbrnq7kZ40LVvizA/BAmBYaYyEGU+OdUVi4BCSBnTqISFuXT6G8QlnnlCa7W7Yw3DLltXEsM6RgPh3w5G2g9JXFDSm4+bOPognJxyMcy8cp4RnBsi7zdFBbx8WryxGtHKwFzM7GGzAPRl4pt/sh8ErL3Z/v6N1ZPfCv5uedovfG4OWI/l0ag39KBLTqi2N7P9ntr/0oKioq6t+v+DUUFRUV9a+WGPhzOzSirTKSb8U5j6ZfsGVRqtUhRt0heTBF+WGHHn2sluCbd7rOqcK1PmKx1OJrxLK3EIOBj2mT7n0nVdmnCxnxnPQNQB7DhwDdu6ChrKrDFoNMZ/8vXbsOijMxC/FaTIdGm7BRySTmMeOy9rQq3JXSDC8Vd+xYy1tEJuseu5XV4TyOI6EXdiKrzN6tZJ+c1Lcj2/3Ic2jQmH5vDkkSVVIir1rtGpJkbApsYxtwn7dCtt5LQreVGthWT4eWat2geFZNA173AWRc2ko5thsHqsZt2BK82LMNtDd4zpmDfDEcWGFXSMaSnIsOW0BpAoJnbwNObjGAnbyVy2pe52Tg3HxOkq123ha7mhWOdam+pegZzasjeld9lyPalrhs0jqAHNEPYAYqD3dEvxV4v25+vv5Xozzw6YjW/0J/i8bgyvpZR/toUtErzX400Oiv/SgqKirqX674TRQVFRX1b9dXW4IN1RpCofXVfGdgZ9Xi2UsVY1XzenlnUtFlh67OMGtK8GruPZdOfNQHA59iWkPacmU7LYH6ykba7c0n1M/FDDMEmv2spZ7XgnwjrosKas/1EoExHrZi+zHl64mamh4MrN0Ss6DN9rxYWpcDX3HydOiTVCsuRaPvNtQXDtIMpR3jlA46damUwhIL1kYY8i0Hi9qmyPY2JIkcLrK93asKLkUb8JYaDS044+XLp31y1pFPVOJ16R5KQIaE3H1EMKdCQ7AW2uWWpB7bBdrrALcdPIywayymK0YLInDBxKmhvGu7xlJWEziOUEGVemqAay3TDt3TPhVJRnHSbn9IwRoXLNn2Un0bHLpKbitstNZ9cUd0z+lV9X3niH40vXscdP9Q4NV7fb5g8A8d0b9dY/BPiYCOioqK+m4VABwVFRX1H5Ts0HlDoNdScy7SocezifdcVFu8JfjVDq3BtpcUvLUE61v+tpi7HXprCRbTQsebzCzCtG/dPjnpYW9m+jQmOfEIKwQaa6hXE6lfArlUE5Jsqf5NvnGq7VRWgZrjQhdI0Lm5TXqdVCyt9eOunmh/ZaOv0qFXPrObmY9TRNeykWDxl3MEsO0MKZUxWkMNpQMtsmJgz1uGNFtgGu5X3pV+1Dl16cq7mtODkGvFBCnfmQvbnwjg1QNarDYuoV5twIW6MR5zGMAbj/fina5sDIY3mPo2062YBZa73NpIRiZY+qxgdGZ3xGgTepfDGenQuRXc8kOPANBtLpv9wNOBa4DWkoJ5YzUn6TIIUEjH29YuNkYUGxOzdnmW+Wx4La5vQuxyRCvkeW3c2dg9/AMaO8dBPdl4PTDaBV6x8WPjwuC3jmjsjF7kX7ox+KdEQEdFRUV9twoAjoqKivpv6ostwYYCUFzb1dP4+KmBDGQqQ5FpIl0/boQR5yxhLnlPLcE7A5/nxwyM1GLAhhtZN3uzxzi9JFfRnXt4E+wKgebLsH3mJC91146c02zxna8RaEmCXvq24rUMtJYmqe3eKkzLLnA6b42+lHwxDXgOLuJjhUzSK4V+Zs0KVjIW1OXZ66sjuwqafbAQTrUPSZpzpJA1BYX16ij2i0aw2bFv1N32W8cTCcIbsdAuWbh4kqsh+Q5YnYXxGGDUihMvTNH4keT7VHyalPp1ZX3HXW3o1E0EWmB/k817Y2CycbGF4yN56DL98YpdLPThfHUFzxws6d5q973Csej9VjiWwfNlXc6XI9qfBTzIdtwc0crWWh85fXjUV7zzKvzX/TYbTGuQPnwB8xfMz/vGlaS1O6J/wcbgcBRHRUVF/YWKX51RUVFR/2Upx/WLdmgXrEZ+/emyQ6/kofVjzcLZ7dCrJfjBwKI1MLAipzeK40QlhBhpy7I3X5OTthZfwIYBM/VMjCle62XGkuucZNcVHO3CslAO+rK/hvIwR+POjlns1iFk7pOEvEWWBOtxWaUstgS/QV+dOzeEYCNUuST6uy+ReQHzitGqHAfl+jByu7wDWZothipPJrSfQqfFpKd8m1RMTkWrcrou3BXjmn1elDK0NQpoYIqviH0MvuUDIGcYC83TmD0PoZ1M0SsSXBCuKCy0VR908uI/+DceMXSAJWYPVZjDYSXulGdTYggUMqINk6XV91KSZv0i2Av4p0TuKzycXc0u1NvLs5PtuE8G3h3RUH0zGHj/iEr1vQ24nmws+z0mJG0b98nAOjIsEu2Wi+6WZtD5hwLvvnFh8Or4XUlav2xj8M+KgI6Kior6bhUAHBUVFfUfF9zKnzLwOQ3PJdUVgLTEz/VT5xDpbHuXb6N0mZvD12z0VUvwewY+ktuXN3szFL9tpG2tTqlgm1QfDCzDtrYXenbXWg3IfELS1uL7PPV9ejCAZyqoWhToJWeIjS+twmiRPebg4tNXCIjtPg7XddHMKcAF3b+lwmC8LwMEmKvaetVTuvRhjVly3XhcVm2NIDY2h5Y7Q7BOGZQzxXl7J+pxGYmZm4W10SQMv+0gOU/3tZzPmU2/Qjun7nQqWRrG+IFuXqPig/inicT2FySFHUN9vHjHGkOwOJuqMxq6rxlOg5OZj0ZaxC3ivWest3GfJgOzhRvdu2yiXm7zc0nBfBAAlG3roUmX6gspuN3Ilvd2c0TTJt1fen3lBdAznasLl0fACOJ+a+vt/ZKCX+OgPxJ4/y9HdMcA5F+oMfgnRkBHRUVFfasKAI6Kior67+sv2KF3e+f10+ye5NeWYODibofeWoLfMrDPfS23bts6J/TsLKqX1HyLuRJ7DE5I6qS4ywjN0GdO3zmkfC49DayZ5qm383LG77HjK+5TB9S9bRXGpGFFN23NzGsYkrqje2nZNsBbTH/vxvaST41g6XzO6/5gGbWt5e1DkuBb5qDgQfH2VGo0uqExJMmQUxuFu1gY8LatOOuGOG9g1fJUa6wxxgsdGIa01oB5SMi6Ypfs6P7PwtSukZWkpXwp4O5RNNoXbxM/AHZsMOpx0Kg8beTDk7FadQZOmraV8dwBInUa0rdhKOCJLil4utnLvDOXFLyNC36ovmJjYHy/ke1yRK+3cm18CLxK0n5ArBqDoZFu1WoAACAASURBVJPfc7BWt/BbgfeL3cLJw8baL9IYbB/8AOCoqKiov1ABwFFRUVG/Sn2FgWV4LsO/kV9Nkrsdmgzs9tGNgTu9nG8YmMN5djlLzltJoM7Ak3OgfvYsdp3H8Jc8Yq5WklOD6FsunhWsD7VwptXvuv8E0vF2KOFrxlAmlyv9m3/OEBpxVeWR1AV4z2U5mc9ptxa84Qi4CLtZYxGOEx2QHfHL2hHsOkcWayV6ymD/hcA9dWPpw0qQ3lFcWGin21OjoQ/njE7jY0vMotDtPdLzkm3Bgvn1cgi/VIYhQnbmTtc5SZhtzLhFHMikaGt/ZICGWrY96/NViPc9qSsYWxokXPRotyLt13uAS2HyNBzgqwd46eHM07qeMjQ0SldNBt6kYFtWlcB73Rls7ErMejqiU37r5F8YvAPzq+rLxx/5zUa+xXvm1ltH9HnvFl7AvLPxL9IYHBHQUVFRUX+tAoCjoqKifqEaP5qQdC7Dc3Y79OPLtzcMJ7dD77qZ20cTiWlOfBXH6nv8KwPb7h71vNmJ4batHZ5kOIHPa4qv/jGZU/LywhVauLvzXrvET+UYrYQqte6yj7TtkrJwejBqGJFOaV4tY728VXhs+F3Ab8o0FqQpKWrBG+KqwGlaOsOu2OPqmi1d0EbXPmaplUV69h+4jnOFMLvpxmCk3CW/X4ZqnQ6TnvLaqFuVWt016vOcAdrzBq7cb7QJt6xs7cGW15Mjf312kd3WI7eWpaVXgiW5t0kExnTlRjd3PeyfeHZAzsxgYiAonoyMi4G9VVhSMP3P+iDppi3hlwvIPhl4uwmtQw/XVay3tY+mzyTE6nH7QC5HtHPp8ClWHzmir8HCE2LRFbypvsqIxgDhcSG3di6MLnvtFv6KI/p8sPHo/21jcERAR0VFRf21CgCOioqK+rXqKy3BY87jNRjAkNX+ZNe9JfjGwKcDdGXqMhhYYHier1rWmCnTnlc8GVgY0MfYPcnL4WxM5BHNZM5dlYUROhVHvjkiadwnJ40ZMY0ryuyJnYro6vIFAh9UiXtfq7L1YP3Zbb3+EzxNyD63dhvSqwiqY4vRqpzJhCG6HPADj/Rc6lJTAZMGUI0zb0sWC3mv78bGxX60ofjpE62qz/5dku+5mb2VKXXC0mz7NDT+ZrUfn0rtImwjrSoB/yAUQ/4jwZaiN0Lvix2+HGip5fAhnwkMnM50Ncv+3d1A7l3B+YCYDXfzZGBonPSMyzZ/zJFINC77IOVNCr7mQq/3eknBZbqRh8eGK+H5ok1XqN88rHnM95IjGsz/wsZvyVZx0PgfyAax+pyn8aHA+0U21kY9sfpPGoMjAjoqKirqr1X89oyKior65eqLLcGIoSpQUCWLPRgYSVKpOQP3m1A8AIuQJC8GpqXYncAvDDy6S8qLNkW2PAvirLoONebqEelU0AvKhOir8bKfy0brSyT4ye4L4sqL9+bpj2OFJ69u4TG7fLGdg390fAQj97yczAhzqkxm5txazeNdi8HM2uNQAPIeo1UZBA3X8dGXTL0PCgYiylosatrmAwH86J2+DNXSRTkBCrbejZat3E1dXV52zbw2zAFiV7AvmK8vjGPWPKTc0loAbNIEb7wZfF/wzEJveuKA5QN01inCexpWwmWioZeCsL29GIMkP7PA+BjAWjKw/qmRSErq8rlQU76mslpTbvtsqrGNRNqlYPAuHnfAon/FTY1L9d012+WI/oiNb/Isb86j11eqL/3z77qF7098lsB7icl3Nt4P+982BuN/ZAHAUVFRUX+p4rdnVFRU1C9aX20JButxZkx96rerJVi0cPteLql2gWt1Wltf6xcAaCowJNlyMfBCr8GYK5w91avLl/JrK7BJqxV5uyjowK1klz1n2rPmGGlCkvH2FQKNObTNkA4HNLyZA4FPGYMp7a6ZSQ7Smk9rQEY9eE0VToZ7dqTm/Nm7oRIQHalRLzFawN1t7K1z5uWdJv71LmTS3N0rCzolzSU25F64a6yLzCd2ya7D2p2DbNnbSo324be6yYq/1oJzNobU+2I3HMFXXMAKG8OAYr47ZMgCR3jBww/6l8fyMONBR8FP7eeKpAbuV0dNZ+CCbO2KVlc7V/YJSYeHP380EkmX9ioFr5FIeCCyYLW56ouA6HaR7agzOHp26r4Nx/Lt+Q0bS/VVC/H1GR6Nbet/UfU9NzbekXs/QkPE97/XGBwR0FFRUVF/uQKAo6Kion7d+kpLsLy+EOgIP4+v3Wj+JJt5htAGkDh645RgIQwZGPrcCwMvW7IsrLsnWftoGjDSnsp+YrdJa2bPOhRys9o2OWkBKnOtNSFJZL7O7RzF5KddTzbQMaIF+E0GlkrcOZ3IcA5RXqWtgcAgk57En4YQx1F606wgxGs9enexztnLWocHeikIeo1ZQrIYbwLYmFHJNUEI5QCoeoVC2ZtYbElZ+LRwFwOZD5ff+4wEkwt6maLd+cxrvGng2wBhnAv25kMB0eiJLZcp2oVrjkQ60BrsRmhlX3dCO842vGUasv2C3gRllfB4iIGvC59S8NuRSAuM16V54NYWjiXVV47oC1Yb5fL8oeq7wyrex+KTgXeybTCR/39ke/vAfzoq6e0RliP632kMjgjoqKioqL9cAcBRUVFRv3TZV+qCJN4ftASLgSEC34cA66drttBzPNKJr9IGC9ihOrVquM7RihHO1qjp5lbjE/AYv3zrq7/vYHiQGvCpXkeGlsh2UOxWb4cCTmsa8Jz6q8G2p8uDB7y7jS9QLhb0Ujqc2zUDWSxt8L2Co6Qkaz0SZiGN7gle56EuX6NY4xRJxB52pWFIc52O8ORA59J5EDh+mZJ1zmFI5wyCBqgefUVbX05g9dbO27WGJD2WVzVGON+eMsjc+9DAwazEVPeZjzGtzs3YKDfMT5KtGmrtMgKgYdVHQ6MruKONXKboxcBYpyum3gMMnzyCwA7/qIiBX6XgaX5eI5H2HO81Egndsmtk9ByJpEcG10f2NfZ5fnTesjHAuD25dKm+f41s16ikh+q7dv7IEf0vNAZHBHRUVFTUX64A4KioqKhfvf6PluDk5JBeIp2ROfUBA3fosdkZWM7SDoy0L+43TpgMDGhJDjxr1pHvUCGEGl5tZ0VnLJOIb2wpyVa6rhB3qal6LcOreBCdWKdg9LFs3+dUvx3bYPbNosG1MwTJzcl8UiKGDgw2LkXgV7tW4y3QM0bLT5szMqEV+Lw9hgCzH8BX++muG2vMEgKm7f/CTsxTE3LGPKXj1p6KntouM/N1ZLBm2tO5nKBmBvKYfbOwFtdiAAkRFQOBh89DGuwj3nqYE4RSkDbWAbs62ow166jzoyMjAN7clYNVi3qAZW5XrzigtEw79OE5WBjKdTc/+52vxaXgmyO6SNBeXcGDPA4wzrDx7ybnnutnbHxeiWJ8UARE34cq/VD13en67c4fHUEPJj5n43+0MRjv0PaZiYqKior6egUAR0VFRf0eJQb+XPYx0DD+XHNrlvH4nISs9KbniGB+d+8ZXa8XA9PDCQbepgQvBl5tw4M/vaBungVUeS7Hsw/vUdOpoGUNEHacpsWXUu8Fk/iHHbu1JQKf5MNSr/GzPoTJFty8vXaHTO2gMDBv/WXSMtZTcmLac+XsIiEopG/cCJ/AhEVmhFpJmCVc+/qxKIYtHwUzfu/6J7peYfyuaZdAMbe3FvvPhW1k40ffLBKqR1KE8p7ONZpzNY4gNs5Ad3CX3Y+eCucqaR7SSeMxCDknWsLZRHxkDimiGNyy+rT92YImA+8MnId7pPMhuRjn7F09uo+RSGvA8v6sYewjkTTbShcIh3JRV/AejiWyfbayN/WTP8OxHIPbDVbRh9vwAfgoxeqRg/UJ2T6O8JHqu8zP+xFkfn44on9uY/DBEdY/5VBRUVFR360CgKOioqJ+m/p6S3BJ9c0ApDklmANyXrzQ+Lbe2Ap7Y+DGruBXBl4ZXBJd90P5WVJ2Vs5XXtSSjudmnVnsfjwONdbkpFIu8lRcViV9geooc/L/J+KC1jouyFSr8BrOJJ35pIk01yKQc0Kbrleol6VeY4SIx7eIacra+hHs0SPd9M8xXdA93bzBMwt60bsKc4JoZl69ypRzixzLnoxl183XY7ovjkA27mOpx1ieITZvb+Vsq1OETM0f3bb9XHZro3c0A2cgKAO29NlKhZow2nT5KMH+nu3y+OQCadV5Ph9R1FZDR/F1Q6bquydjneRwz/16CcfS5xAO9vlmrZFIu5C7s/E1Q5g7v7IxHuV0TAb+yBH9VrPdnxZ9ovoKpH/oiD7/+cbgiICOioqK+ssVv0CjoqKifqfSl+nP82/EqCA4oNszFkuziDwV6f5dHK7a1jzLamNg+jn7wuknAx8gk2U8Vil869KLk88xguDZ2fqLGKplbebZK/6NhOKthxlmacnMk1EH+3g7rbvCaezQzsXTwPXBRtNle86IwurYs62u3Qx47lD88o3idCGL5XwgcCljj5imzmmH52LLYxjSqflGiROoNkTU9oRpP0UjhR3pybGtZR+SBFXU0DUL6bP8yRKBqZFDFx/OxrgDMwAMJNaKfNFyessmDXcytXc8TTiAtOepRwrscCbs6x2nwxnSLBqtOYRZ46YBxPWKxZIdGv3O9CLvF76ZwC/IHyv367hNQj4/CMeCos1xwfCWbzvLBP7WEf3UjZmD9SYc61PV9/8l2/0IP2Tjn9gYHBHQUVFRUX+nAoCjoqKifrP6ekswtDSwxZuWYBDZ8cYpfTKDqmimkehiqrLihI8YGJrbq9qcm1jXFVqqr52jixdDXt/kx3t7tqdG57KMyisui0JiM5hcIrBeQud2vqYNkaIBfI1YWxAltW6gy6HJBwLvScsIOk5XjJZ+tNh4DQQGr3Lu0YqDqsXjr6+o5GMyf6l2g/EgYM4Exi08rsAw5Vq5C5qIDmyzfxr2Nhe6vdGXg3w1SEkLxhryHP5MXVQzgQGlHPikkUh6OZ4LtOHp0HwvFPwlIz0YWFnQfI9wKG7BQZi4JvLkgKNbDtZTCr7U3bPU9pSCGY7VGW2Fj8re06uRSE9HNBT3xwTghcEPy4Mc0TcxmWdcOVgfpVjtO/811Xc/Qp3jtXmxP6ExOAA4Kioq6u9UAHBUVFTUb1kd1JA+t0OLZCA5vmNg+pf7M1N3/hg5wiv6qOSLgWdL8IOBgR+Yc3s7DmTGJQ9Wd+eKXQ29lI3s+q2qAZgTg5J2cMKJKBfuOu3CVE5Bpjt3isBSa90mzWlDEFZLU8esHR/auOHu5UwG9iv0WIFSfup7jNYtYpqTk6APkz/TGoZUIRqLjZUd7VxKRFS+lHD6iomydzFfgVUSS+1mlnKpu4CrPhuDW1OztM4Iu3u+OpAlAndwX1rp0ANx19BPcbc5Ihh2dN4WXLkhma06M8GbvLs+OeTSUmVqHnAXg4FT8huu93QUScGPHCy76jWraX0kcHU5r9Tu/ePos4LbtXPnQw19RF9x9+GIxsOCTAfBuLXvSvV9OKL78Byst7D6RbIVSH+FjcdPbQyOCOioqKiov1MBwFFRUVG/a8Gy+iMGBt0pzejlG/w5m3U7u0xfGRgKXPIG3SXMgj5mS7A8zIuBkaJ03gBGjb7oHZXaWYpY+mLX5L2mjkKD3bAd7Zq7h/kSirNHQ+1xWScbTR1TpwiMPmH0mbrQalQoEViaKtBFDCZnsjK2cjaYhPp9esr0Oc+Zc5Wsupzejru5Qx++T+jV+CX701J2t/ZkY/dUHxeol4b75slNs28WA4QNVlva96QBu6jLd0/Mhop7pBW4BVrLDHy2u64JSdXld7QWwx49FJct57YunHOcGzd7dreE8UPRZcUZWO/1lQWtRyQ0MAP8+jMHi9Oe6lMKPqcF4M7GA2/8nJP0wsYIiO7bR3Q6olcsth+i12cL8XREv4ZjyS3/ker7FTZW0/hHqu8rG2tnOaL/cmPwD8PwoqKioqI+qQDgqKioqN+49DXa+OtzO3RGz6h/fX+ItGrW7cz7fTCwujEfDDxf5S3BDwZGM/B9XJAYuPUZ9XwgdvnCyHFK8LyagZn5LMHz6iydo4NFyGocfWjaimV+pD37AirikRlM5SRttwtjgKcA6yIzY5ZBqiPtAc523swYLbmadwyTPowlvejGdD9j7tGePIwxUYcbnm8HgS/70nt953GUdKGjq75s692bpWEtPtpKzIJ635uak3Ht58zoYhuzM/DhNmmYnzNs0orRhpae+KYy9FkMXPARU5tu9yxoJWYZ/ldOzyqe9f2/9s7EsHEdCZsMhaEwFITCTBAKQkEofOhuEAQBUJJtzXieWPUfu0tdlOSxXPr6yLOgz8XP9mr0UbA1AOcoeDd/+5GLWrLd6K6tREpmW9ctDxuA5Z6Ds9z4JKtBKgsGbhzXvkTiasKzpb5Wz3y83c9S3/pgXRH9vcbgx+0PAADwGAQYAOB/z4stwZNK7BpdL7qzFOmuta0dd+6tLjWpqrcdwvlWe0vwMcsq5tlIsvJ3bmuYk+Mk09BBvDI5+TTpSrwrZDmKxyJi25dTpkbnR9EYWULg87gsM+1cS2yLl6pUNj2E6KitxpWcMFnkJCXNez7plmhZaB5YHaejqTUeY7RkgnM4Gb6NwrJ649NxzYf1y4lTy6uN0WpqgOVsnW8OWifztifVcR/Q5RfRXZuMZb4kjb5TMtJF69NFSZeqplpalL0qnqbfuhZKvhHQSmkJidO7L/Xhev2g8bB0UC95f5I8eAj5nnVMdHFgqdNeZbSYfEEw7bOgt3YtcH6C+jXESXd1ard0X3dunF44UXS533Ac1hp4+Rk+627R4NOV0yl7qSN4McjVhVAn666vXHvp49S3duOtKp9u3Lgs2S7fYbnVvdgYzAhoAICfwO9QAIBPIL6wISkXPC/TYA+weo7uQBo4sGXIct/JD9ZDO6P208rA4fVwYEn7Qh47XAqJS1AsRcmL5rfxtK03N+KmP/4l4tTOz22v8rV7F/3eynhqG5XcOGceg5VrXZcyKys/fU1l5SJtoZx1FFbOJ7UbWc5wzX3CYonrvqE36Z+N0Zok/q3dOJ+5Vl/LcTfXKa4+6yDdyMtcTz+2RcE6hOwYlSxDraZ1qwqzpew5vcA25kr1MofAQRcIeVd2AlsILNOh5PuKVZq+58VC4PzckyMvroi0PDsrkw6y/lcceJb9xXs8nhuAa/9P5yy7lBZpAF5FnKX4WbcoHbOg5e3bJ3VbFNyuBU4/RG4uKX1560136x8YefWCnKFURIdWd+WNcOdaBv0BklcgnoZjpffW+cGMaBeXpuQh7t3CTW78OMht/x09vHIT8BaRbiqiHzcGeynmXoYXAQDAKyDAAAAfwistwbYDSUf7ngcLWSWt0503fvDHd04C5Q94ucLhwHtLsNPq2XweTmJPGflbCV4dFKdHsEXE7UPY9mAzXV0vlCRWtNnLfCObqSTx4X5fUtDrp9o5LQSWnth5rVt266cv4ieJ29F4afOoktzJSt55X0SkS5tso48ln1ZBLYYTF9uoZM9OarwXmZ19LP61xl19RjKEzC/1ca2p9jYQy7YZ5yLnRd6avDXX7jld0Wm8qYXipfA7PZyVNMtTWLQpWm+bG2K1uvh0JuL5u0hrP7BNzxLb16pjqUK3JmpVUHNgZ1XAJQfWqVcimT43AEtdgPUPi33Ox0KpGM36BrOgRbR9joLPw7Fykfy6Ng3AFgU3Ghz1a5RGd+30+opoicTDlN6hUzYb5em1g6N1jbAUn59l9cGE515Wr1Lf/sqNG7/SGMwIaACAH4IAAwB8DuUP6KctwclZbDLw6c/xfXnS4mU17OCGSS+T+IXQSPIaJMuSGLByYFG1ZBNuMMtKNCP9kb+eapi3vSFZim+ds0Loo+I6aHTs9ofQplNNGpO2ulMeO2+mrF6qvefmiaisJiubdJnyUQ0uN1z3gc82LHrLJcFWJu0nl3cpbSEvBN4Tzvy1wjxpfujqLbhyaovqkN35fnxxUU3SHdXCsg5IqqnTeyPduDq7K72AeR+S99bMvOi3D8FnUZcUUZNnsV8Nge20xaxCHtCVn6Cap+3XtTJpWS41uTw6K6mrGGWs7T0dkVhYZ1zFPBtr7waPEgXnBmC7t221KDhvjdpD7EV+zk5Ztz2c9qS3UXDs9iTZBcGGn1UpuvynVj43Rc5a3762g6P3iuj0f+vc+FGQGwe6axFx48YPUt/eje3n6rEb143BzT+QjRHQAAA/BgEGAPg0XnFgpyuCrBS5GYtlldLen4K7ckOZ9pskQjesNn/xNw4sA7G22fLeclfFgS3EE48Kg0e3TFFac+tpz3pX67TmQmi9z3Sp9PqGo+O0bP3VsDo2Q7lMVpOkOl1hW0fEoohun6QVguXP266O6bhOzM4rhWy2c2khljbpxdsy5GMLrji8TnQKk42kLvtyk8Eu6u1xnxpt639L767OT5ZR0FLPvax1Q6+8Zctc5FZicBnTnHcy1UOe0qNIo+/+HYRk6bqZOXfeLumM1dV3412SOc5lrFkuKY8aQadzTM9ltUFdyfndrCG/DgabY/k+wiZL16q/7QH7YBa0bgAeDMdyne7aeznZMuRBA/BpcLRasDzNrvLZ6vPrSVpbqYh+uci5tO82V7Yf2Rev3HcLF5He9m+yrKmh+Ydsbhz2PmoAAPgqCDAAwAfytCXYVNbN2W3KBGPDSozDcunAyX6SrSWRaPI0WbGrDiw32oulLazLDaim0Pu8KzPwdFf1/VvHqRXE2raeY/BvTFInpif34H2Jas2mZMLWXpNcJi2LHJ6HcmmauKzz4hfXTHWW2HParF1WBliNBjtLBfXeXWzHxSinxULgEsCWFmLJrYNPirVZH7I6sy5AlsDWBlMVN7aO31yHrOKd3sbJ68Cv0lC9yf+2KVZ2c0sU5UwmeRHqLlaZCjZNViZtz1AKof1q31BIMn32Uvmp0Ag3zw+L2YG9fmEgX3ysuxW7xWxcFkdNmhXrQ1gUnF5dOefl9AI+mgWd3qyzG8u72V25NAA3uisn3HcF29Q338WtQTyzHYd+XbfcLwl7vch5q3p9G929ioirZ5ALOup/yPY/JR+WKd9xAwCAL4IAAwB8JlKGLOtbn43FSiazjFuCxVbc2krIfsNkg16KjE+lrdYSnITNlEoUbg8kpYpWU8FcSm030snSMle5zgbNgd1muW5OkpvTnnRwsRQ0+2PY1S51McQcQqofShXxeqiXWu4cJt3row2fzVRn6ae9ntQlAmkF3vvLZQmnTYESo95HT5vWJg+UWuJ5Ku24FvlqlL0URQxRjif5zKOq7RFNvMPecrwflFprWSeV09ryCqTXQsS46q+OujZJPNDlhl59IyT2Lxms1F3LffnyCqRzK6XF+VkvufJ51dLvMjRLJmlrFCzfStQNwKPUV/YbW4dzo7sXs6Dlm5HzZOxt110rzK6vXGZBn/pmZbmwFNv3umtjn7+UzX7pyuX1LFce6u7Tmc/2D7ns/rUR0EM3BgCAV0CAAQA+lldbgpcg033XpR0IZKOhnS/Jak2ulF7muss3XxRC7cDH3SZVW/b9tMWBV61K9ae4L+ZJzzpmaZr6/UzWLWx2lBynZHHygJpRSkHyXiRsZdIifvtwKQnQ1B5NV8oG3fKswyrBaboXGbNcvXrlca3a+Xitdq2tFwJvuxtLr+uyNHGo9NzOwULg47hUMcv/O5VtbyGHwPtB0c502lrhbC240hgsEi32rgH3ceZmvBLVirlqd/SSR0bbay2dxekJrYfxeptTvX8tErVvXDbxqkZLObQOJDO/zSoY51LgnR/X1gJXHdH2HKUBuIuCYz0LunrBQ9RRahdRcLMnyZxZvvzwZ88M/sEG4Nd7euV7B50c3hy/moPV33O6Wt/W+wpRR9wl+61HQKPBAADfAAEGAPhw7E/kB4NzxCacjsWap1YSVDQWCzD3BLJg/bR+lsixbNzJN9Ry6NqBcy4X97lKIWYHtolZyU/8Wpt2cWCZ/7RPPz4/r71f13LISjttFZD1xMrN94lWMv1Yh0slkbAsV0ZPbTJoyo5nc7O0ddGW2CiW1ZyVd3L/83kJkzxN3bHbNDbL8XR1LZBulE8iX9fanVUyNyoos7WWfVHwvs5XYkzVXesftpcoT9LSjDqZqkXuXsNlOWLmOk3mt9l4dU9ymOfjSOXAVi1sDizF23vlsywino6xz6a71gBcD8GyrwDqcuj8MloUXDX6bnkW9KXulgFd5bhU6Z8naW37++fXwbwr+aoltgYrTdRhbLCv1y0/uHIdBX+PIOvH5F/xsPs37kOzaAwGAHgFBBgA4PN5qSXY7w4cl+FYrFxD27UES6Y4r3nxb32RhaK7Axc9iOpb4sCWIqoDyw6hpAoWMZakNO4PmBRXe1DPpbC5SDvObauwXZzsMMuYzsrKXbvpmS4uCbCX0u98qcikVeeWKc1iucGmcM1+qffZ5l1Ki4ynOnprq4uS7UlR99l1tah76eceJ60VMT4ft0XBUmu9n096zdKDSWXyftBeUh9kkrPEuTLRKofPq5Py7xikOFlmPru8Oth6X1VuxUjFeMNSAkmZmKUObO+UObDcalms8lkHS0VrrrZHtyFYeeyzVT5vazrPJvVN78uSTt6+Yqji7qso+NDd+ofNdLcvyLdty6PhWOll0SnV1c+MNgDLy9JVRMvL2vyAXetuv6PowZWbL25eJO5zsPTLmtn+ywPFldXWNAYDALwAAgwAcAvs7+nHG0RlZrIT9ygWVC5KlrGnjK2BWDaocWnnwFFrgNV6vPT7Vu2g2hJcO7BtRRIvmqsyWgtjxeBnKYT2rSyJqWqZtLTEnluFN62HtULo2sCTRUhFcHKKKchCXZm8HPIkrTVvBjL3lht67fP0vi7utW3D0mob2iFScpGOnq7Lqrd9NrWUmk+n4053NclrXt2PPrrMj84htsakouKL7jFa8sms0dYl6yvpg0tvjyaZ7plJtwAAIABJREFU8rKXAu8423cT8rz0JvmZRmkhtlW6JTy3V7J1YFnvu6+lXffKZ3Vgq3xOb33+okGrDHQzU3q3fF5xbG/FHgXXXyXYM3Vl9dG50dfrKLImY9e53oMoOFdEr66+E/kSxMsQ7GbSeIi5IvqkpjHPu2rcOI56erdtK7rbp7798t4XibLDScJek95VKRddzYWub05FNADAUxBgAIC78GpLsA0Tcq0DW9grDtwVnZqIHg58vn9xlkXMJetsuVRbgkWgdge2xlppB63WGq2mUTE3sg5CVxvIpH25TYOllSVbS21pv/Te6zJVGU+VZ2V5XwJkiSU1q3Qy6yoWz5QhxrMrs6x1LpTuKJJq42PzrZ1PslExw+p4folWCZPL/djxdN7pSbl5Xwisx9NLJfXGiyublmQrcvrfGrfaRGUZ+BxkI7EdjPuSJJmHXExeJ2YlgQ+z7hN2zszNHtoks/6+I71Pmy5Jts7tw4GDzzOf/SIF87vcSrobZ5sBVufVZd71abxZtTiqrXwe6W5YQ3bjke62UbB+S9J/QSP18Lq0qdHd9DPWjDG345KKd8eTM1/19DZR8Ddowt5VNkuP7/AVxUWDAQAegwADANwLyT+f/AEtDpw0xXK/diyWZGRL3ltzvpV2C/tm51B9w3SF6FsjkoU6shNILw3eVuZu21bPrFptpZGYqy4unk93Yl3KktYucnbJVdqH1rFVll8H2aPkcxIug4X2rT9qOOUmTrNKyWydLM61G+7fDojmWaKbvyYoHa16tqbLsmw2eXV1XFJ0zY2T7MzuOL6mR5DocbKp0XY8vTI52k3mKW+ZFBUnBRdbnvZh2mtuRU53U0ZkaQB7RL4y9lnbWUXs18naniXq1NRXvmuYXZP6yrCrZOk+ljz5mOBdVT6fZj7vQ7BOteIWBXdDsIb7kErF9aDyeW2vbIdtFvRLe5LkuxbR3SYKtnRXfqLP6a79MPgw1t2+yHn7OlHSaV+HvQ+69Pvb1nOhH1znGycGAPDxIMAAALcj6kTZxw5s7bg5Y2zGYmmldK5EbcLetaqUHjnwsobswNXxpCXJfKSTVw25TMyKIbcExzI1Ws1cJvHqJK2yGjf9qS/NnrYidesyQJ/HZdl+oHNlqRZR64gmGyuVj2s5dHI8a/f1fimLf+141uNlNWMvx22+lzz7NZdeu30hsM39kq8P5uzMx/E1SLPurrV23IdjkLW4sRYVL9LpGYp5Jq2dvEsHnXmmJofiaWGWG5gD6xEbX5xOxqTRXC6qG6bnIu2yukA4G++8WlYs07yT/Is/Hg5sUfCynLLcrLvn1Uel8jkXcpc3vXpZ2qLlLgrWdzA8qHyulyodx/tu4agF4f58J/bixEG6m34sh7r7vTHO+ykcYe/jeuZX7sp6g2n6BQD4EggwAMAdKX+IP7hOGYtVt4mWiw4HPsdQ5sCL7ypXd9uc9X6bbuGgG31kGdO6mXjkity9JThPjXYyMtrqfkWKllOBsTW+Slusn7bzU8uWqxuVZANSdc4WEdeDlMtFopezdBHbGOpS6Wrl0OmEZR1x1XtcyqRlNNQSbcWRSaCVQ6fnMmlunHto9+vLiCi/HIHquqWXIgmnm1fLbPNz1O28SUqLGNtLZKlsHkZlZ6JSpxumcvpqFeDWRmsjpko/sI2/knbZRTf9Vg4sLhmCeWBxYHtQW/+bH3RfcSzO/LjyufqRGEfB+lXCMAq+rnweOHP6gZEfm6UTbE3C+0Zf+aHsZkHniHh79M/kKfZvrQ573zWumWpnAIBvgAADANyX5y3Bqi4S7HVjsSRjc/vw564lWIZMubVeIFQuEotZB1OjRc2W1Ry49KnaTWQNz5KnRosdaU2ymZs/O7CU7zotanWnWVNRJ3VpZe88OVHI8+uQ1+faIK7mhGXC1CRjqOuLLKqVe1piHZLnymR149UfRlrKpNMTnN2hgiFJpZZb50na+4lZ5Ctdynveu+2BqvRn74XQx9sUZ7v+IaibpNY21+rkwFrkbA6sm6ryK2ltxun1sWVI+Wp6Putmk8DSq3oUOZ8qn88zn68qn2ddRNzrbhn0Vb/uV5XPg7XA23YMxxpFwe1a4L0iutFdawBuu4V/UORcdhc9nsH+Q9BgAIAvgQADANya+HRDkjnwPJks1d2SInXJsNZ2i8+2e+MjB7Ybnh3Y8klxYBfrDT3b3hJsDizeq823dnelJThnudYrG+d6kta2J73pyDxNyzpYU+R1I5EMkTrPNLIbhlnGULebb1Qkk6rWA6hN6uQM5mBhb308ObM+j7mOcJ2s6VU/XKpFQerM07556HQ+aphNyioP6NZTSKst3CXy3XTQsbyq+oTTQwbNSC3cDmvQ2uz0n8H5ySqfiwPHMgx5cbV7p5+HSYdg1d2/RXeLkNc/S33lc4h+OATLe+l+bnJjud0XK583y5mb3FUrovu1wDI06wdjnF8cZ/V24r4NmIpoAIDHIMAAAHfnpZZgLRyVcHU0FstGQzflprHs6e0ceGu6hc94GYo8cGBrEJVuYe2wLWmt9sJm87FNvEkC7ba2//awNS+lxUmA5yjhZxMR62Dn5N/zsItYHHiZ5q09W6kQ1t2/TXQsUW16gsmw5qmp+9UdTN4KpPOjmzEuPpeOl2dsDbdurZ3TjpeJ0I2LuiSz+/cRZq0+BptrZW5cO7BM81rno/t3zWO6bQqUVFZHWdN05NLa/SslxJOrxV6On7t/NxsQ7duD8WoIVomCL1LcXmvHa4Htjexb0P3a70Oyiuh2H9IXifvuIgt7nXO/mMS+PkmrYdp56+kAAPyL8JsOAACOKsoH8VFpCZ7PS2KijsUSc9tjxtNFWind7C7Kd3jtwEkMk6mKO8dTnaoU66rLiayuq61NsovK1GiT1SJ1dUScH1RdPbmiD6fUscizlDSv7QlLmfQkM66TLzVypYKmm4SbUluNjvsnKMeTaFuJb3Pcr3VonE9MU215eas64aSXk9Za1/G7tv6qc055ZFTUDmHJcmep4s4tytr9u0k/dfRBM8+9+zfb5h7AilFHPVjl1VKvHPJdlfOxcujmoOnu4gaVz271eUbX+QcmfzFxfvHTOYbZ9Slu8Fr53JVPXzmzSPbSvk3f4xfD3rdTey8ODAAfD7/mAAAg89SBfTjGYplBHbf1sh4p92F26mjdwjJL64x1C+c8sInsNMlzi28dWONoWf+jDixbkUpVtgiOlBbLI+4OLOuLtt2ItEw6/YmfFDG3v8bQRMRyw2UWG2xczhJTHfssW3bOFdTypJ3vS8Fl9ZFseBp8NTC7MDxuBefNcR+DDbKup45ZaCwznPfgNNuplyrlpJdlbrNsOY66KzhppxqmNXUHXQQlN4vBun8lNpfXxtuDi9bqhOg8t6y0Ipvuns9fHij9YHSNvhIFjyqfH0XBXYp7VD73ke9VRLwOKp/1ufrtW5Swt0jvB7Td9saLAwPAZ8PvOAAAOHi1JXg0Fktsylnz68ADpTNXfamdMGRF1P1GnF1tlrl14E3zXpsanRSwrE2y07O8N1lYMsa6aNmmRkfnJttgNC25btm0uYqIJRdN10p3q22xTRexc+nWU63Nm6XHUpjt+hhc7Gzx/cRseU0sG+8M0OL0RiNXCWQPrd2KcOr2XZsaXXRXNg/pYGfTeKkXLg6c/C/I+WQH9j47sEa+1v0rs7D2J1G6f3vdVZNue3etHNo2S9VPqoy8fmUIVvqvwyFYW9mT1DX0DiufN+8HB79IH/Z+Up8tAgwAd4PfcQAAcOIVB5axzHt1bxPPLjbheWrz3qhzmCUH7mpQzYHdPHDgrexV6hxYRjrN3hy4VvGoi4VtF5Ht5DnW2EQ55cUWC8+Sn9rIqyYilq3Cc+4xLpXVdgcWEcsA4ej8+YTV1SU9zg5ca/O8ScHzwGklNpfrn4NNqTZfpzAt7ajnJYi3L0djrTlweunEgctiYdVdU1Np51W9T89c+oF9eZpBcl53cuDNNj9tctAmRJd3wwZBO/umYNf43P0b5ubrAKlwjtm9m/lVNgq7r3Aer/9Nqt5HwZvM6erv5EEUvH2d+Md2F/07mOg29c8Nv3d2AAB/Cn61AQBAy9OW4HosVu/AyVvyaOiuitWtuYi6d2CJXh86cLLgfk1rnhq9nBx421uCiwMfrcLeO3WqqHmwlFLvN7RlSxKiarVzOpkkoFJf3XURLzI7eJbhUueUUp6CCzIxWyuTawe2Qugg86fXknbKAzn9UmCemgx8DRIO+/NkKf2mwMLbkI9r5CuvjG5rXnUPr+muTMlSqywpq80wCz75u8uRrxYbyzBknf1VIl9pRE4Hg6vfjRwFx6mZaqYSPu7+He89Ss+4d+Oy/vc82zlXPg+19qrRd/Xbd6nDXvsa6JPC3prit30DMN4LAJ8Nv+MAAGCMOfCDubJlLJZMCo6nPT12kTlw7br1RcnE6nurHbgesmVYM+3AgfUO5zk41zqwlcsmITSVtRuq3jhzJ1lGpDH2Mds5RuvvzR3I82r3eWiY91mPXbAW4joiLuZsT6F2QisCT04rZ1Ppnzqwjr12rq4QFqsMYtchBCfSmzUyxByVBwnqcxCqVdO6nUjqvPNxc2C5Tx2pJecj1eLSF12nu/JuaFF37cBb1f2ruf3+8ljkGyffaa1UOHdLoSUKDq5x4+1odV56fR1GwfI6j74c2dcXnw5+FZNeJ7XxwkeGvQ216CLAAHA3+B0HAACXvFQO7RZz4PQfpww26LbbqQ0MNzVnsbhJdgs3SmMTs9Idthtr9hy1d+Btnxq9rMt8XlNktc21A5ctNbYhyRxYFc0dD2N6HKKtcTq6iNeY+2xNdBefBybHUBzPyrnTa5K0U8Jwt5Sa4fzUpPY61PonN3H5Ed16RMdqj3M5XrpqvT2c3n+JgmXMs8WtIZTjNrJ7VV/fcgIvRcVSslw5sH6J4KRWXIvCyzcecmKbtnRXDryVyLeb4OU3nyd4vebGzXcE5UEfRMFDZ96+RSlzSJSw1/YYfWrqaxS/bXS3/y8AAB8Jv+MAAOARJR+7vIK1BO9jsco8KrvIjKhfg2QX2QioxnVLKOr9YGq0FOyOHNjaensH3nT8lcwmTiaYBMdVmbYtFp4Xc+BynzHmxcLpqNxnWGs9Drt8Jq+UVbS7W+ZJWuqfefKzqnLR5q00/c45Oj6cNshsrTKh6tDpLW9RaiqHZQ723h5chkuJ0lazmu344cCljVm6gZ0Mu1qPCme7f2mpHkW+k7xua/0jYMcvtbYb6x23k9vX79r6lShYv2wI23fpx1k1YW+x4g8Y7/yUXoDp/gWAO8AvOAAAeMKLLcE2FmvZTnZaLpJy2uUUEe91sHJ84MBz1CbYwdRokdbQzcTaw8OhA8cQLQeelimcW0+31RYLh+Y+j4h40Q1PpZBbyqTzfCmZX727Ze5iVZfT2cOhaH922mWRVNmp0+Ys94hAc2V4jpS3on+imutUNLtEx/qqHhmyDZeSu6j3Felx6fINU+3AGvnqFxbrucv3YeQ7revZavPxfuj3MPLd9jLsQfdv5fyn4xYF/yyPjV/fXZRu4pybtf78Jw/9j/M4AcaBAeBT4bcbAAC8RJKBx8lY8FntBmOxkuLIWCXdoHseC5wukg7YPiLeS2H92rqu9f1KybMf2FG6laha58CbTGiOIkGxlWqbNiwO3GTLuiHJ9Pg0kFmW9ITSRVyffHoNpD9WCq11YNVy9PpaF2t24HU6uauqbHbg3XXXfdZxiL48ummhlDcvs/drvRIpH5+nkwPrcaeLierzjLr0SGqem7LnPfJtKt/zcXXg5riNQxvoq65r3s61A1aGPYyC80Sxc8/5t+X3adj7FO/94/r//zvDmuer/w4A8DHwqw0AAF7FmiSftASvEglK0fFZQcM+MStd0oyA9mt10Tkitot6B96snHh1QwdOjzV0YB1tlXxvbcdl5ee2JNGV0VPnml6LiCU5PstfEqrSRXwS+Bit9Vd6crW593BgjYilslprvGt3NcWVumhX6fQejSbfrg28lDdbCXel3yF5cZjnvK9of2VsMpbEw+v5i4YYXVf2bKc0x2UY+crx85Xl/reQK8D7yHdtFyBv5+j7dDz4n4Su6T1M1lqHvQ9GuL14hx9cET1cgHR1BQCAj4FfbQAA8AVKsHZ5hS0PkVr93HimmJ4I2zRYgxS0lnia0q1i8PVFJdg8GnF3ktdlw+z8RFzXD0RXx/xuyRKHDhxEB0PfkmphZu/Am25IyhFxna/uEfEqPcenhtjkeH5yjQPbg+dWXr82d2U67eelSaGtvNlLifZScnWLWNPBXFC9P67F5otE09NWvX0S7Wqc3pc9l8i3cWBpXh4dt1r3vtm72aVczifKRK7L3vLX6cPe906xijoK7kELwP+Up36LAAPAR8KvNgAA+DKPW4K3fadu0MU7p3LoTZbniiap67ZFszoaWuckdxI1r7bCp28JtlsNHXixHPiQxmh/09swLRHkdgaS7cD1w7FMQRwrNNXaeTpxMi/d0Ftf5CWJFge2AVd1r61pc9LRekS2uqskul5Gb52O52rw5MCzO5Uxa3nz4HjwMvTLtet20/sy6cBnWwVcMNftC35l6ZFuCm6Oe3nVff/Nw7Acequ7f8/H65FpX8Kktw57/3S/7udFwU/zXgQYAD4SfrUBAMB3iE83JNkQqW4s1rbrcb8GySZm6fbZqVkFrJmqOLBMb+pagtNFtnO48StVylhkMjlSya5tztZAdHVnr2TOXWip1uV7B940pg6TvB7tzid13d6BLSI+HLjkt3s5dOPA2142HKalOeegDyHXb3bzynalIKfUH/ey+1dHkp1e/CU6WYM7++aLCZe0eW7HL+couCuHlu84RoOdc+T7XePd7ySHvWV30U/u7XuP/jHDsR4bL/YLAJ8Kv90AAOCbPG0JzkOkkrWG1oGTHqvjOblC1xKczNBGQLd6rBOzhmOx0kUSwI4cuOyV9d7X+5yOi0bpsXNrW/BsF2lELAXGnR7LHqFzqLvtEfEgB961OR2vHVhvEsxp+7TZu6V3YJuYJXH71D60dErPbVG3lT3n2uxuZ++yDqLaR8eXthx62yPf6N4gqD8fZ/VeovKLJ/AuhorLGiQA+Hj4BQcAAN+nZGJXSpCHSO1jsU4FzLbrdZL+03YNUtiHaXVrkEQaxYIGLcHZZkfWmkuIp8l1Q4ll7tQyuInaexyI7l53PXRgiYid61th8yjmzoFNm7MDt0574cDphM/9wOUhZqfDpbqlRDqPa7yKeXw8tI695XR68PLq/ch3BekEmuNtl/VXsN1Fs+Kkdflzao8BAOAXQYABAOCnvNgSLItnu7FYm1YvJ9GM7jRYq/T9Julq9NjMedgSnCcwj1qCxQ/npV3js18kHaqjyV7ZgbuLpGDY+bK1qLnJsqx9mfSmo5h7B7YSblkTtbR5rLju5Oq50PkhtBzay5cAnc9rxXVzwhJoJ0vtd+3GfSDz0r34wyVG15OuVh3AvP2Mfy3sBQCAzwMBBgCAN/C0JVh7T3NLcNPfa64lrb9Lq2dSHa0RcT8WK7o8Fqsth96HafX1ui4xjfU4RJk13d9ETs8mVPfF1VpInIctd89VapVHy42Dd/29mevKLqgui87l0FOb62pEvMjxzmlXHWnVP4QcHz0RWao0mlwlD31ep1QuyMffQfrJsbC3SC9hLwAA/DkQYAAAeA8vtQRrcjvHqd3rq3osS33Wbg2Sl+VD4y3BOjFr2BJs5tzIm1bS+qzHfahr5jzs+y0XdXp8OHA/f8uFdu9u6b+t9gPXt1lcO0o6P53VZddt7mrdpH67d9ogC4EH5coXo790Wthg9rWZ+fB4+FkfLGEvAAD8CggwAAC8jS+1BC/b0tzY1LQfi5Ujx25z0rY39162BNtYrF3ebH6vmHPwVxOz1jAo8d32jlm5iff9Resy0GPr1x2WSUdd0SR57/nerBy67y7OrjuPNgnF6EbLkOWE3chpL0Z/WXmzPET31YAcX7pq868T//ruol+HQVMAAP8a/PIFAIA3Yw7sO1EslJbgponXBiYPx2KJPLlrPb5uCZbZW3vNc23mVxOzngyUjt2iI7voWo99DJcRcWxnX21VOfTYzy3X7bR5jaMxV/pcpBy6d1qpbm7Lp+02MpnsTeXN+10eYa/VyX/GFOWnDC2XVUMAAL8Lv3kBAOD9vNgS7GUmVVI6X19i/undoOY5aDw72BK8B8vJgZt7K8FyerjWNx5OzLpS0BwsD0dDX+txtsr+3vxgq/B2XXStLb7tWqOt1ubOdUO8vKvg3NVXANvPaMLe9P/frbPXftiaH7nHPgwAAH8Bfu0CAMAfQSJb5cE1DmuNbaib7Ld3YLkkaLfwIruFm/t70BJswbLIRjfI6mpi1hYvW4LDRd/vttv7UGgl0x5WUEdtFR72JPtxbXN6lH40dL5oNLZKHt2PWoKtxXcaPJHvUWrgE/U3IN77J1+IfBZFaxFgAIB/DX7tAgDAn+J5S3DczIFlSU8/FksvSg48XIMkZdKdHluwvK7dkK38ONM48LxqCd7ylqZHCW1/UbjcnPTliPhiBlV+Ea5cN4xPbLjBKD/NH/DiOKvyw/Dxfb81CDAAwL8Gv3YBAODP4r1/PB1a9uDMg5ZgYc01z40DH8O01kF6XDYk1fcmIaQWXV/VPF+2BJs3jlQzjhLa7VpoN9Vjf+HAl/cWg5Rq95O0zM+/4rqmzW8ZZ2W7iyzsfX13UbrVPaPg4f8cHgEAgD8Kv3YBAOCP88qGpOFYLLto0zVIMgL6jCW33nVbgmM253ST0hKcBFjqseMe6l7NfxrmvWatwxlXlh6PCpivhPYqvN22y4hY66QHI6Af3+SN23r3u3zD7qJbRcEIMADAvwa/dgEA4G9Q3OnBNXI5dGibePXGzoR2qMdSDr12ChqiLRa2e5PtQnuIKnmyOXAzTnl7qMcXNc951vSVnV5VHV+Et1lcR0Z9Vdt83KQv4b6exf0i6bX33te7i96S394kCkaAAQD+Nfi1CwAAf4/HLcF2DXFgP57zbKuA+5bgctHAWledmLWJ/Z50S815XD98vQYp63E/2cu2NL1Jj/MzGlZQJ63vVvhutVG/iT7sffvuojtEwQgwAMC/Br92AQDgr/I8+tsHWbmt63qN+4ak0N3cFgi7dktw6Rael6ldTZzubhnPuDqs9UG3cKfH8tAX9/Yg7w1uEN5uFlNf5b1uXA4d0v/5Ac3uou9VOH/jQZ+Ux/+fQYABAP41+LULAAB/m5dagrXmOYnmsOY5ObAPoxhWO3XbsVgqk/M8ubXdEnyEuqOa56sR0IeCjrqFc7D8FaG9yntzTH2R90Z3XU/+FUrYW3YXvT3sfeUEPjIKfmy82C8AwN+H37wAAPALPN2QVFqCZUvw2VrrbuHeWs2c+4hYHqtqCT7dKOyhbr+k14/1+MGGJA113fDestAOG6GvI+JLPf4B9vo75/5m2Pv0lH79HN7OUHGnnb9/PgAAwC9fAAD4NZ62BEdtxx1aa1mDFKI/3WSveW62BMsDJb+yW3URcR4BfaWgF6Fu7hYeTp/6htBe5L32ZMfHv0i9u6htigYAALgBCDAAAPwmUrb8WkvwwIG15jmJbi+0kiYuU7q0pMdH4Labc1Nc/TjUfSC0Vnc9GFgVHgrt38p737K7CAAA4DNAgAEA4Jd5pSW4lEP3LcHRrcOa59wtvMpYLHuI/qK2uNpC3dGMK+FbepzLob+Y96bbDY6/jJUTv313EQAAwP8dBBgAAH6f5y3B23blwKXmebwl2PJet7QriHdrHZqzhbox+P6iRzOuRnr8Ujn0myDsBQAAeAwCDADwI64m2bDv5BuYA7fLimpUTZPoDsqhQ655bh1Y79dPk7sOdcfp8VNrvZpx9dW898HzfYH0kqRXrA57H72AAAAA94a/xgAAfsorrsvQ1xd5pSU4uiU58HCQlcS2y+Ci5IROtwQP1dQaiXtzDmsY1zxvLzTx9tGrXTQcAf116rD3V3YX/S5Tx2+fEQAA/D/gAwMA4Kc8FWD7n/yN/iJiuI8HFEvoKZHn2IHXQUvwqjt/TE2Ht5LlSctlzfOVOV+FumLOV3nvDzTVpLcOe29b4cy/JgAA+B58fgAAvIH6z/Gh/fbH4QGvtATbjCcXxjXPUti8HllrMmorDLZbyXbf5hH3Tt0rqV79Rc3z9U4jaSR+RypbXo2ERb5MtOJfEwAAfA8+PwAA3sPTmJc/2b9KUtbH06GTtQY3jx1Y/bPsOioCbBel48NbxaDWOgx1gzjwoFtYlhVflEP/gAfjrJ4Pzf50+KcEAADfho8QAID3gAD/CZ7KnjhwWJc4WgWczDSuyXXTFWaZ9hxOF/l12Pcrme2yJD0ehrriwN4NT/TnAlx2F1nY+6AO/KWh2Z8L3b8AAPBt+NgAAHgbj/8W5y/171Gy0AfXEQcOo2HO25YcOFnrPE2NK5o5j0NdEVE15+j7i9btwpy/y7d3Fz0fGPahXHUZAAAAPIXPDACAt4EA/zmeBp7JSJPNujjy5BjTiz9o7lXVdWHU97vpXiU3rnmW6/9s+lQJe4v0fs9jbx4FF/jHBQAAL8IHBgDA20CA/yjJGJ+WQycHXuIgnhXP3C76fqWHV3YLX93h0IG/Rx/2vkVczYFvGAUX+McFAAAvwgcGAMDbQID/NK+0BPuoLcGbLwfNnEtL8NCBrwqbo3T3tge/es5/YXfRzSdj8Y8LAABehA8MAIC3gQD/BV4p+jUHTk6b/6f3pYU4O3C3Bqnc6oe1zc15pse13UXvCnufPuIdyqGfbt4GAAC4gg8MAIC3gQD/NZ6bXkwum5t7kwA7d4xutprnoQO/JexNj/VHw97HpEe8w2SsB5u3AQAAHsBnBgAA/C95OgPZCpul9Vec1DUXyXToOFpo9K0zKbuL/gX5LEXXv3safxrWIAEAwDfgYwMAAP6vvNL4GsI6uentNvhJuGAAAAAFK0lEQVTt3UUAAADwiyDAAADwP+aVxlfn3OIfrRF+/bHesrsIAAAAfgsEGAAA/veYA3vvh5cuy/KTeJawFwAA4GNAgAEA4BN40BL8jcHI6fpJp+uw98quAQAA4H8EAgwAAB9CslbnXO/ASWJfFOA+7P34lUIAAAC3AgEGAIDPYdgS/HhQcJmZTIUzAADAx4MAAwDAp+G9L9OhbVJ0f51/bXfRL8I+IQAAuA982gEAwAdSNiSFEJLfloOMs2qovRcHBgCAj4ePOgAA+ExKxmsBr4W9pr6/fWr/Cr3x4sAAAPDZ8DkHAACfTFJfm4xlYW9Jhn/7vP4JEGAAALgbfM4BAMC9GA7KuicIMAAA3A0+5wAA4I7gwBsCDAAA94PPOQAAuCnJfpn//PQIAADAJ8HnHAAA3JebtwQjwAAAcDf4nAMAgFtz55ZgBBgAAO4Gn3MAAAC5Jdh7/9sn8ldBgAEA4G7wOQcAACDcsyW4Nl7sFwAAPh4+6gAAADLJgZ3y2yfyV5l2fvtEAAAA/jh82gEAABxYS3AI4bdPBAAAAN4PAgwAAAAAAAC3AAEGAAAAAACAW4AAAwAAAAAAwC1AgAEAAAAAAOAWIMAAAAAAAABwCxBgAAAAAAAAuAUIMAAAwBuYOn77jAAAAKCFj2cAAIA3gPECAAD8+/BpDQAA8AYQYAAAgH8fPq0BAAB+CvYLAADwv4APbAAAgJ9C9y8AAMD/Aj6kAQAAfkojvTgwAADAvwmf0AAAAO8HBwYAAPgH4eMZAADg/SDAAAAA/yB8PAMAAHyNV9p9EWAAAIB/ED6eAQAAfkqvuwgwAADAPwgfzwAAAG+gNl7sFwAA4N+ET2gAAID3wBokAACAfxw+pAEAAAAAAOAWIMAAAAAAAABwCxBgAAAAAAAAuAUIMAAAAAAAANwCBBgAAAAAAABuAQIMAAAAAAAAtwABBgAAAAAAgFuAAAMAAAAAAMAtQIABAAAAAADgFiDAAAAAAAAAcAsQYAAAAAAAALgFCDAAAAAAAADcAgQYAAAAAAAAbgECDAAAAAAAALcAAQYAAAAAAIBbgAADAAAAAADALUCAAQAAAAAA4BYgwAAAAAAAAHALEGAAAAAAAAC4BQgwAAAAAAAA3AIEGAAAAAAAAG4BAgwAAAAAAAC3AAEGAAAAAACAW4AAAwAAAAAAwC1AgAEAAAAAAOAWIMAAAAAAAABwCxBgAAAAAAAAuAUIMAAAAAAAANwCBBgAAAAAAABuAQIMAAAAAAAAtwABBgAAAAAAgFuAAAMAAAAAAMAtQIABAAAAAADgFiDAAAAAAAAAcAsQYAAAAAAAALgFCDAAAAAAAADcAgQYAAAAAAAAbgECDAAAAAAAALcAAQYAAAAAAIBbgAADAAAAAADALUCAAQAAAAAA4BYgwAAAAAAAAHALEGAAAAAAAAC4BQgwAAAAAAAA3AIEGAAAAAAAAG4BAgwAAAAAAAC3AAEGAAAAAACAW4AAAwAAAAAAwC1AgAEAAAAAAOAWIMAAAAAAAABwCxBgAAAAAAAAuAUIMAAAAAAAANwCBBgAAAAAAABuAQIMAAAAAAAAtwABBgAAAAAAgFuAAAMAAAAAAMAtQIABAAAAAADgFiDAAAAAAAAAcAsQYAAAAAAAALgFCDAAAAAAAADcAgQYAAAAAAAAbgECDAAAAAAAALcAAQYAAAAAAIBbgAADAAAAAADALUCAAQAAAAAA4BYgwAAAAAAAAHALEGAAAAAAAAC4BQgwAAAAAAAA3IL/ABl8zOcLt+N6AAAAAElFTkSuQmCC","width":1400,"height":719.1257,"sphereVerts":{"vb":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1],[0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0]],"it":[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270],[17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288],[18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271]],"primitivetype":"triangle","material":null,"normals":null,"texcoords":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1]]}});
  rgl.prefix = "";
</script>
	<p id="debug">
	You must enable Javascript to view this page properly.</p>
    </div>
    
	<br>Drag mouse to rotate model. Use mouse wheel or middle button
	to zoom it.
	<hr>
	<br>
	Object written from rgl 0.98.1 by writeWebGL.
	</body>
	</html>
  
  {:/}
