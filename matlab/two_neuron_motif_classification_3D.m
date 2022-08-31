function motif_class = two_neuron_motif_classification_3D(x1, x2, y1, y2, t1, t2);
    % Two neurons are the same (two distinct neurons involved)
    % Assume t1 <= t2

         vec = [0,t1,t2];
         n_distinct_times = length(unique(vec));
         
         if n_distinct_times == 1
            % Two neurons are the same, all times are the same
            motif_class= 4;
        elseif n_distinct_times == 2
            if (0 == t1 && 0 == x1 && 0 == y1) || (t1 == t2 && x1 == x2 && y1 == y2) || (0 == t2 && 0 == x2 && 0 ==y2)
                % Base and first nodes are same; or first and second
                motif_class= 6;
            elseif (t1 == 0) || (t2 < 0)
                % Synchrony is first
                motif_class=7;
            elseif (t2 == 0) || (t1 > 0)
                % Synchrony is second
                motif_class= 8;
            end
        elseif n_distinct_times == 3
            if (x1 == 0) && (y1 == 0)
                if (0 < t1)
                    motif_class=11;
                elseif (t1 < 0)
                    if t2 < 0
                        motif_class=10;
                    elseif t2 > 0
                        motif_class=11;
                    else
                        error("Shouldn't be here")
                    end
                else
                    error("Shouldn't be here")
                end
            elseif (x2 == 0) && (y2 == 0)
                if (t2 < 0)
                    motif_class=9;
                elseif (t2 > 0)
                    if t1 > 0 
                        motif_class=10;
                    elseif t1 < 0
                        motif_class=9;
                    else
                        error("Shouldn't be here")
                    end
                end
            elseif (x1 == x2) && (y1 == y2)
                if (0 < t1)
                    motif_class= 9;
                elseif (t1 < 0) &&  (0 < t2)
                    motif_class=10;
                elseif (t2 < 0)
                    motif_class=11;
                else
                    error("Shouldn't be here")
                end
            else
                disp(x1)
                disp(x2)
                disp(y1)
                disp(y2)
                error("Shouldn't get here.")
            end
        else
            error("Invalid number of same times: $n_distinct_times")
         end

end


