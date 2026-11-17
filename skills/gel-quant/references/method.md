# Gel quantification: method and reasoning

## Band signal is the area, not the peak height

The amount of stained DNA in a band is proportional to its total signal (the sum of pixel values
over the band, after background), as long as the signal is in the detector's linear range. Peak
height depends on how sharp the band is: longer strands, edge lanes and overloaded lanes give
wider, flatter bands with the same amount. So area is the measure. Bio-Rad Image Lab calls it
band volume; ImageJ's gel analyzer measures the area under each lane-profile peak.

## Lane profile

For each lane, take the median across a strip a few pixels narrower than the lane, row by row.
The median ignores dust specks and scratches; for bands that are even across the lane it is
proportional to the full 2D volume.

## Background

Rolling minimum, then rolling mean, over a window wider than any band (41 px at ~2 px/row band
sigma). This removes the smear under the bands and slow gradients without eating band signal.
Check that the window is at least 4-5 band widths.

## Peak fitting

One Gaussian per expected band, all fitted together, centres constrained to +-2 px of the expected
position after the lane shift, widths 1.5-4.5 px. Joint fitting is what separates close bands
(e.g. a 90/80 nt top/bottom doublet). Hand-drawn boundaries split shared signal arbitrarily.

Every visible band must be in the band list, including ones not used in metrics, or its signal
leaks into a neighbour.

## Gel smile

Lanes run at slightly different speeds. Each lane gets a shift = (anchor position in this lane) -
(anchor reference position). Use the brightest point of a band that is present and strong in the
lane, or the centroid of a doublet when the single band is absent (e.g. lanes before the last
ligation).

## Signal to molecules

SYBR Gold / SYBR Green signal on single-stranded DNA scales roughly with length, so
molar ~ area / length. This assumption has not been tested on our gels. Test it with a known
mixture (e.g. equal pmol of the 160 nt and 90 nt strands) on the same gel.

## Bounds from clipping

If a band in the numerator is clipped, its area is too low, so the ratio is a lower bound. If a
clipped band is only in the denominator, the ratio is an upper bound. If both, the value is
unreliable.

For standards: each clipped standard gives an area/pmol that is too low. The largest ratio among
them is the tightest lower bound on the true area/pmol, so sample pmol = area / that value is an
upper bound. Valid only for unclipped sample bands.

## Size calibration

Mobility on denaturing PAGE is roughly linear in log10(length) over a limited range:
y = a + b log10(nt). Fit it from markers or known bands in the same gel. Extrapolation beyond the
calibrated range, and modified strands (biotin, spacers, dU) shift positions; treat those
assignments as inferred.
