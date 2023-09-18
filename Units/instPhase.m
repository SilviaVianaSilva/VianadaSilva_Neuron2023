function [signal_phase signal_amp] = instPhase(signal)
H = hilbert(signal);
signal_phase = angle(H);
signal_amp = abs(H);
end