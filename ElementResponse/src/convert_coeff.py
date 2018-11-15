#!/usr/bin/env python

# Script to convert an ASCII beam model coefficient file to a .cc file for
# inclusion in the library. Whenever the beam model coefficients file are
# updated, this script can be run to automatically update the corresponding .cc
# files.
#
# $Id$

import sys
import time
import re

def flat_index(shape, index):
    """
    Compute the flat index of the element with the provided (N-dimensional)
    index of an N-dimensional array with the provided shape. Row-major order
    (or "C"-order) is assumed in the computation of the flattened index. The
    index is range checked against the provided shape.
    """

    assert len(shape) == len(index)

    if len(shape) == 0:
        return 0

    assert index[0] < shape[0]
    flat = index[0]

    for i in range(1, len(shape)):
        assert index[i] < shape[i]
        flat *= shape[i]
        flat += index[i]

    return flat

def regex(name, type, signed = True):
    """
    Return a regular expression to match a (possibly signed) int or float, using
    the named group syntax. The matching group will by assigned the provided
    name.
    """

    expr = None
    if type == "int":
        expr = "\d+"
    elif type == "float":
        expr = "(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"

    assert expr, "unknown type: %s" % type

    if signed:
        expr = "[-+]?" + expr

    return "(?P<%s>%s)" % (name, expr)

def main(args):
    print "converting %s -> %s (variable name: %s)" % (args[0], args[1], args[2])

    HEADER, COEFF = range(2)
    state = HEADER

    shape = None
    coeff = None
    freqAvg = None
    freqRange = None

    line_no = 0
    count = 0

    fin = file(args[0], "r")
    for line in fin:
        # Remove leading and trailing whitespace.
        line = line.strip()

        # Skip empty lines and comments.
        if len(line) == 0 or line[0] == "#":
            line_no += 1
            continue

        # Parse header information.
        if state == HEADER:
            match = re.match("^d\s+%s\s+k\s+%s\s+pwrT\s+%s\s+pwrF\s+%s\s+"
                "freqAvg\s+%s\s+freqRange\s+%s$" % (regex("d", "int", False),
                regex("k", "int", False), regex("pwrT", "int", False),
                regex("pwrF", "int", False), regex("freqAvg", "float", False),
                regex("freqRange", "float", False)), line)
            assert match, "unable to parse header: \"%s\"" % line

            shape = (int(match.group("k")), int(match.group("pwrT")),
                int(match.group("pwrF")), int(match.group("d")))
            assert shape[3] == 2, "unsupported array shape, expected d == 2"

            size = reduce(lambda x, y: x * y, shape)
            print "coefficient array shape:", shape, "(%d total)" % size

            freqAvg = match.group("freqAvg")
            freqRange = match.group("freqRange")

            coeff = [None for x in range(size)]
            state = COEFF

        # Parse coefficients.
        elif state == COEFF:
            match = re.match("^%s\s+%s\s+%s\s+%s\s+%s\s+%s$"
                % (regex("d", "int", False), regex("k", "int", False),
                regex("pwrT", "int", False), regex("pwrF", "int", False),
                regex("re", "float"), regex("im", "float")), line)
            assert match, "unable to parse line #%d" % line_no

            d = int(match.group("d"))
            k = int(match.group("k"))
            pwrT = int(match.group("pwrT"))
            pwrF = int(match.group("pwrF"))

            index = flat_index(shape, (k, pwrT, pwrF, d))
            coeff[index] = "std::complex<double>(%s, %s)" % (match.group("re"), match.group("im"))

            count += 1

        line_no += 1

    fin.close()
    assert not (coeff is None) and all(coeff)

    # Write the output.
    fout = file(args[1], "w")

    print >> fout, "// Beam model coefficients converted by convert_coeff.py."
    print >> fout, "// Conversion performed on %s UTC using: " % time.strftime("%Y/%m/%d/%H:%M:%S", time.gmtime())
    print >> fout, "//     convert_coeff.py %s %s %s" % (args[0], args[1], args[2])
    print >> fout
    print >> fout, "#include <complex>"
    print >> fout
    print >> fout, "// Center frequency, and frequency range for which the beam model coefficients"
    print >> fout, "// are valid. The beam model is parameterized in terms of a normalized"
    print >> fout, "// frequency f' in the range [-1.0, 1.0]. The appropriate conversion is:"
    print >> fout, "//"
    print >> fout, "//     f' = (f - center) / range"
    print >> fout, "//"
    print >> fout, "const double %s_freq_center = %s;" % (args[2], freqAvg)
    print >> fout, "const double %s_freq_range = %s;" % (args[2], freqRange)
    print >> fout
    print >> fout, "// Shape of the coefficient array: %dx%dx%dx2 (the size of the last dimension is" % (shape[0], shape[1], shape[2])
    print >> fout, "// implied, and always equal to 2)."
    print >> fout, "//"
    print >> fout, "const unsigned int %s_coeff_shape[3] = {%d, %d, %d};" % (args[2], shape[0], shape[1], shape[2])
    print >> fout
    print >> fout, "// The array of coefficients in row-major order (\"C\"-order)."
    print >> fout, "//"
    print >> fout, "const std::complex<double> %s_coeff[%d] = {" % (args[2], len(coeff))

    i = 0
    while i < len(coeff):
        if i + 2 < len(coeff):
            print >> fout, "    %s, %s," % (coeff[i], coeff[i + 1])
            i += 2
        elif i + 2 == len(coeff):
            print >> fout, "    %s, %s" % (coeff[i], coeff[i + 1])
            i += 2
        else:
            print >> fout, "    %s" % coeff[i]
            i += 1

    print >> fout, "};"

    fout.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "convert a beam model coefficient (.coeff) file to a C++ (.cc) file."
        print "usage: convert_coeff.py <input-file> <output-file> <variable-name>"
        sys.exit(1)

    main(sys.argv[1:])
