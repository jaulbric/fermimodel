from xml.dom import minidom
import numpy as np

import Tools

def PLspec(roi, radLim, maxRad, varValue, var, sig, nO, f, i, p, dist, TS, vi, fixAll):
    """PowerLaw spectrum"""
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    spec = xmldoc_out.createElement('spectrum')
    comments = []

    fscale = np.floor(np.log10(f)).astype('int')
    
    spec.setAttribute('type', 'PowerLaw')
    
    comment = xmldoc_out.createComment("Source is {0} degrees away from ROI center".format(dist))
    comments.append(comment)
    if dist > roi[2] or fixAll: #if beyond ROI, shouldn't attempt to fit parameters
        if fixAll:
            comment = xmldoc_out.createComment("Source parameters were held fixed in 3FGL analysis, free at your own discretion")
            comments.append(comment)
        else:
            comment = xmldoc_out.createComment("Source is outside ROI, all parameters should remain fixed")
            comments.append(comment)
        spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
        free = False
    elif dist > radLim:
        if dist > maxRad:
            comment = xmldoc_out.createComment("Source is outside specified radius limit of {0}".format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        elif vi < varValue or not var:
            comment = xmldoc_out.createComment("Source is outside specified radius limit of {0}".format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment("Source is outside specified radius limit of {0} but variability index {1:.2f} is greater than {2:.2f} and varFree is set to True". format(radLim, vi, varValue))
            comments.append(comment)
            free = True
            spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
    elif TS < sig:
        if vi < varValue or not var:
            comment = xmldoc_out.createComment("Source significance {0:.1f} is less than specified minimum for a free source of {1}".format(TS, sig))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment("Source significance {0:.1f} is less than specified minimum for a free source of {1} but variability index {2:.2f} is greater than {3:.2f} and varFree is set to True".format(TS, sig, vi, varValue))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = True
        spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
    else:
        spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        free = True
        if nO:
            spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
        else:
            spec.appendChild(Tools.parameter_element(1, "Index", 10.0, 0.0, -1.0, i))
    spec.appendChild(Tools.parameter_element(0, "Scale", 5.e5, 30.0, 1.0, p))
    return spec, free, comments

def PL2spec(roi, radLim, maxRad, varValue, var, sig, nO, F, i, dist, TS, vi):
    """PowerLaw spectrum, but no flux value from catalog"""
    fscale = np.floor(np.log10(F)).astype('int')

    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    comments = []

    spec = xmldoc_out.createElement('spectrum')
    spec.setAttribute('type', 'PowerLaw2')

    comment = xmldoc_out.createComment('Source is {0} degrees away from ROI center'.format(dist))
    comments.append(comment)
    if dist > roi[2]: #if beyond ROI, shouldn't attempt to fit parameters
        comment = xmldoc_out.createComment('Source is outside ROI, all parameters should remain fixed')
        comments.append(comment)
        spec.appendChild(Tools.parameter_element(0, "Integral", 1.e4, 1.e-4, 10.**fscale, F/10.**fscale))
        spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
        free = False
    elif dist > radLim:
        if dist > maxRad:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0}'.format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Integral", 1.e4, 1.e-4, 10.**fscale, F/10.**fscale))
            free = False
        elif (vi < varValue) or (not var):
            comment = xmldoc_out.createComment("Source is outside specified radius limit of {0}".format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Integral", 1.e4, 1.e-4, 10.**fscale, F/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment("Source is outside specified radius limit of {0} but variability index {1:.2f} is greater than {2:.2f} and varFree is set to True".format(radLim, vi, varValue))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(1, "Integral", 1.e4, 1.e-4, 10.**fscale, F/10.**fscale))
            free = True
        spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
    elif TS < sig:
        if vi < varValue or not var:
            comment = xmldoc_out.createComment("Source significance {0:.1f} is less than specified minimum for a free source of {1}".format(TS, sig))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Integral", 1.e4, 1.e-4, 10.**fscale, F/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment("Source significance {0:.1f} is less than specified minimum for a free source of {1} but variability index {2:.2f} is greater than {3:.2f} and varFree is set to True".format(TS, sig, vi, varValue))
            comments.append(comment)
            spec.appendChild(1, "Integral", 1.e4, 1.e-4, 10.**fscale, F/10.**fscale)
            free = True
        spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
    else:
        free = True
        spec.appendChild(Tools.parameter_element(1, "Integral", 1.e4, 1.e-4, 10.**fscale, F/10.**fscale))
        if nO:
            spec.appendChild(Tools.parameter_element(0, "Index", 10.0, 0.0, -1.0, i))
        else:
            spec.appendChild(Tools.parameter_element(1, "Index", 10.0, 0.0, -1.0, i))
    spec.appendChild(Tools.parameter_element(0, "LowerLimit", 5.e5, 30.0, 1.0, 1.e2))
    spec.appendChild(Tools.parameter_element(0, "UpperLimit", 5.e5, 30.0, 1.0, 1.e5))
    return spec, free, comments

def COspec(roi, radLim, maxRad, varValue, var, sig, nO, f, i, p, a, ei, dist, TS, vi):
    """Power Law with subexponential cutoff"""
    c = (1./a)**(1./ei)
    f *= np.exp((p/c)**ei)
    fscale = np.floor(np.log10(f)).astype('int')

    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    comments = []
    
    spec = xmldoc_out.createElement('spectrum')
    spec.setAttribute('type', 'PLSuperExpCutoff')
    comment = xmldoc_out.createComment('Source is {0} degrees away from ROI center'.format(dist))
    comments.append(comment)
    i = (i if i >= 0 else 2.)#some pulsars with index1 < 0 assuming standard convention, means rising counts spectrum at low E, odd
    if dist > roi[2]: #if beyond ROI, shouldn't attempt to fit parameters
        comment = xmldoc_out.createComment('Source is outside ROI, all parameters should remain fixed')
        comments.append(comment)
        free = False
        spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
        if c <= 1.e5:
            spec.appendChild(Tools.parameter_element(0, "Cutoff", 1.e5, 1.e1, 1.0, c))
        else:
            spec.appendChild(Tools.parameter_element(0, "Cutoff", 2.*c, 1.e1, 1.0, c))
    elif dist > radLim:
        if dist > maxRad:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0}'.format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        elif vi < varValue or not var:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0}'.format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0} but variability index {1:.2f} is greater than {2:.2f} and varFree is set to True'.format(radLim, vi, varValue))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = True
        spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
        if c <= 1.e5:
            spec.appendChild(Tools.parameter_element(0, "Cutoff", 1.e5, 1.e1, 1.0, c))
        else:
            spec.appendChild(Tools.parameter_element(0, "Cutoff", 2.*c, 1.e1, 1.0, c))
    elif TS < sig:
        if (vi < varValue) or (not var):
            comment = xmldoc_out.createComment('Source significance {0:.1f} is less than specified minimum for a free source of {1}'.format(TS, sig))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment('Source significance {0:.1f} is less than specified minimum for a free source of {1} but variability index {2:.2f} is greater than {3:.2f} and varFree is set to True'.format(TS, sig, vi, varValue))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = True
        spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
        if c <= 1.e5:
            spec.appendChild(Tools.parameter_element(0, "Cutoff", 1.e5, 1.e1, 1.0, c))
        else:
            spec.appendChild(Tools.parameter_element(0, "Cutoff", 2.*c, 1.e1, 1.0, c))
    else:
        free = True
        spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        if nO:
            spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
            if c <= 1e5:
                spec.appendChild(Tools.parameter_element(0, "Cutoff", 1.e5, 1.e1, 1.0, c))
            else:
                spec.appendChild(Tools.parameter_element(0, "Cutoff", 2.*c, 1.e1, 1.0, c))
        else:
            spec.appendChild(Tools.parameter_element(1, "Index1", 10.0, 0.0, -1.0, i))
            if c <= 1e5:
                spec.appendChild(Tools.parameter_element(1, "Cutoff", 1.e5, 1.e1, 1.0, c))
            else:
                spec.appendChild(Tools.parameter_element(0, "Cutoff", 2.*c, 1.e1, 1.0, c))
    spec.appendChild(Tools.parameter_element(0, "Scale", 5.e5, 30.0, 1.0, p))
    spec.appendChild(Tools.parameter_element(0, "Index2", 5.0, 0.0, 1.0, ei))
    return spec, free, comments

def CO2spec(roi, radLim, maxRad, varValue, var, sig, nO, f, i, p, a, ei, dist, TS, vi):
    """Power law with a subexponential cutoff. 4FGL definition"""
    f *= np.exp(a*(p**ei))
    fscale = np.floor(np.log10(f)).astype('int')

    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    spec = xmldoc_out.createElement('spectrum')
    comments = []

    spec.setAttribute('type', "PLSuperExpCutoff2")

    comment = xmldoc_out.createComment('Source is {0} degrees away from ROI center'.format(dist))
    comments.append(comment)
    if dist > roi[2]: #if beyond ROI, shouldn't attempt to fit parameters
        comment = xmldoc_out.createComment('Source is outside ROI, all parameters should remain fixed')
        comments.append(comment)
        free = False
        spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
        spec.appendChild(Tools.parameter_element(0, "Expfactor", 100.0, -1.0, 0.001, a*1000.))
    elif dist > radLim:
        if dist > maxRad:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0}'.format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        elif vi < varValue or not var:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0}'.format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0} but variability index {1:.2f} is greater than {2:.2f} and varFree is set to True'.format(radLim, vi, varValue))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = True
        spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
        spec.appendChild(Tools.parameter_element(0, "Expfactor", 100.0, 0.1, 0.001, a**1000.))
    elif TS < sig:
        if vi < varValue or not var:
            comment = xmldoc_out.createComment('Source significance {0:.1f} is less than specified minimum for a free source of {1}'.format(TS, sig))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment('Source significance {0:.1f} is less than specified minimum for a free source of {1} but variability index {2:.2f} is greater than {3:.2f} and varFree is set to True'.format(TS, sig, vi, varValue))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = True
        spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
        spec.appendChild(Tools.parameter_element(0, "Expfactor", 100.0, 0.1, 0.001, a*1000.))
    else:
        free = True
        spec.appendChild(Tools.parameter_element(1, "Prefactor", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        if nO:
            spec.appendChild(Tools.parameter_element(0, "Index1", 10.0, 0.0, -1.0, i))
            spec.appendChild(Tools.parameter_element(0, "Expfactor", 100.0, 0.1, 0.001, a*1000.))
        else:
            spec.appendChild(Tools.parameter_element(1, "Index1", 10.0, 0.0, -1.0, i))
            spec.appendChild(Tools.parameter_element(1, "Expfactor", 100.0, 0.1, 0.001, a*1000.))
    spec.appendChild(Tools.parameter_element(0, "Scale", 5.e5, 30.0, 1.0, p))
    spec.appendChild(Tools.parameter_element(0, "Index2", 5.0, 0.0, 1.0, ei))
    return spec, free, comments

def LPspec(roi, radLim, maxRad, varValue, var, sig, nO, f, i, p, b, dist, TS, vi):
    """LogParabola spectrum"""
    fscale = np.floor(np.log10(f)).astype('int')

    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    spec = xmldoc_out.createElement('spectrum')
    comments = []

    spec.setAttribute('type', "LogParabola")

    comment = xmldoc_out.createComment('Source is {0} degrees away from ROI center'.format(dist))
    comments.append(comment)
    if dist > roi[2]: #if beyond ROI, shouldn't attempt to fit parameters
        comment = xmldoc_out.createComment('Source is outside ROI, all parameters should remain fixed')
        comments.append(comment)
        spec.appendChild(Tools.parameter_element(0, "norm", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        spec.appendChild(Tools.parameter_element(0, "alpha", 5.0, 0.0, 1.0, i))
        spec.appendChild(Tools.parameter_element(0, "beta", 10.0, 0.0, 1.0, b))
        free = False
    elif dist > radLim:
        if dist > maxRad:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0}'.format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "norm", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        elif vi < varValue or not var:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0}'.format(radLim))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "norm", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment('Source is outside specified radius limit of {0} but variability index {1:.2f} is greater than {2:.2f} and varFree is set to True')
            comments.append(comment)
            spec.appendChild(1, "norm", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale)
            free = True
        spec.appendChild(Tools.parameter_element(0, "alpha", 5.0, 0.0, 1.0, i))
        spec.appendChild(Tools.parameter_element(0, "beta", 10.0, 0.0, 1.0, b))
    elif TS < sig:
        if vi < varValue or not var:
            comment = xmldoc_out.createComment('Source significance {0:.1f} is less than specified minimum for a free source of {1}'.format(TS, sig))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(0, "norm", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = False
        else:
            comment = xmldoc_out.createComment('Source significance {0:.1f} is less than specified minimum for a free source of {1} but variability index {2:.2f} is greater than {3:.2f} and varFree is set to True'.format(TS, sig, vi, varValue))
            comments.append(comment)
            spec.appendChild(Tools.parameter_element(1, "norm", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
            free = True
        spec.appendChild(Tools.parameter_element(0, "alpha", 5.0, 0.0, 1.0, i))
        spec.appendChild(Tools.parameter_element(0, "beta", 10.0, 0.0, 1.0, b))
    else:
        free = True
        spec.appendChild(Tools.parameter_element(1, "norm", 1.e4, 1.e-4, 10.**fscale, f/10.**fscale))
        if nO:
            spec.appendChild(Tools.parameter_element(0, "alpha", 5.0, 0.0, 1.0, i))
            spec.appendChild(Tools.parameter_element(0, "beta", 10.0, 0.0, 1.0, b))
        else:
            spec.appendChild(Tools.parameter_element(1, "alpha", 5.0, 0.0, 1.0, i))
            spec.appendChild(Tools.parameter_element(1, "beta", 10.0, 0.0, 1.0, b))
    spec.appendChild(Tools.parameter_element(0, "Eb", 5.e5, 30.0, 1.0, p))
    return spec, free, comments