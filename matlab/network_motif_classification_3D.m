function motif_class = network_motif_classification_3D(x1, x2, y1, y2, t1, t2)
     
     vec1 = [x1,y1];
     vec2 = [x2,y2];
     vec0 = [0,0];
     
     list_vec = [x1, x2, y1, y2, t1, t2];
    
    if t1 < t2
        x1_new = x1;
        x2_new = x2;
        y1_new = y1;
        y2_new = y2;
        t1_new = t1;
        t2_new = t2;
    else
        x1_new = x2;
        x2_new = x1;
        y1_new = y2;
        y2_new = y1;
        t1_new = t2;
        t2_new = t1;
    end
    
    x1 = x1_new;
    x2 = x2_new;
    y1 = y1_new;
    y2 = y2_new;
    t1 = t1_new;
    t2 = t2_new;
    
    % Assume below that t1 <= t2

    if all(vec1 == vec2) && all(vec1 == vec0)
        % All neurons are the same
        motif_class= one_neuron_motif_classification_3D(x1, x2, y1, y2, t1, t2);

    elseif all(vec1 == vec2) || all(vec1 == vec0) || all(vec2 == vec0)

        motif_class= two_neuron_motif_classification_3D(x1, x2, y1, y2, t1, t2);

    elseif any(vec1 ~= vec2) && any(vec1 ~= vec0) && any(vec2~= vec0)

        motif_class= three_neuron_motif_classification_3D(x1, x2, y1, y2, t1, t2);

    else
        error("Invalid number of same neurons: $n_distinct_neurons");
    end

end

