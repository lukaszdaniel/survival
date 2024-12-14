/*
** Fill in a new time-dependent variable
** 
**  id: subject id in the baseline datata set (integer)
**  baseline_time: start time for each interval in the baseline data
**  new_id, new_time: the id and time point for the new covariate
**  new_covariate_values:  value of the new covariate
**  indx:  starting point for matching in the baseline
**
**  Both data sets are in order of time within id
*/
#include <R_ext/Boolean.h>
#include "survS.h"
#include "survproto.h"

/* First routine, for cumtdc, return a cumulative sum */
SEXP tmerge(SEXP id2, SEXP time1x, SEXP newx2,
            SEXP nid2, SEXP ntime2, SEXP x2) {
    int k; // current index
    int n_baseline, n_new, last_subject_id;
    bool has_cumulative_sum = false;

    int *baseline_id, *new_id;
    double *baseline_time, *new_time, *new_covariate_values;
    double cumulative_sum = 0;
    double *updated_covariate_values;
    SEXP updated_covariates;

    // Get lengths of baseline and new data sets
    n_baseline = LENGTH(id2);  /* baseline data set */
    n_new = LENGTH(nid2); /* new data set */

    // Get pointers to the actual data in the input vectors
    baseline_id = INTEGER(id2);
    new_id = INTEGER(nid2);
    baseline_time = REAL(time1x);
    new_time = REAL(ntime2);
    new_covariate_values = REAL(x2);

    // Duplicate the baseline covariate values to protect input
    PROTECT(updated_covariates = duplicate(newx2));
    updated_covariate_values = REAL(updated_covariates);

    /*
    ** i= index of baseline subject, k= index of addition row
    **  last_subject_id = prior id, id's for baseline are integers starting with 1
    ** has_cumulative_sum: 0 if nothing has yet been accumlated for this subject, 1
    **  it it has
    */
    last_subject_id = -1; // Initialize to an invalid ID
    k = 0; // Index for traversing the new data

    // Main loop over the baseline data
    for (int i = 0; i < n_baseline; i++) {
	// Reset cumulative sum when we encounter a new subject
	if (baseline_id[i] != last_subject_id) {
	    cumulative_sum = 0;
	    last_subject_id = baseline_id[i];
	    has_cumulative_sum = false;
	}

	// Move forward in the new data until the new ID is >= current baseline ID
	while (k < n_new && (new_id[k] < baseline_id[i]))
	{
	    k++;
	}

	// Accumulate values while new ID matches and new time is <= baseline time
	while (k < n_new && (new_id[k] == baseline_id[i]) && (new_time[k] <= baseline_time[i])) {
	    cumulative_sum += new_covariate_values[k];
	    has_cumulative_sum = true;
	    k++;
	}

	// If we have accumulated values, update the baseline covariate
	if (has_cumulative_sum) {
	    if (ISNA(updated_covariate_values[i])) updated_covariate_values[i] = cumulative_sum; // Replace NA with cumulative sum
	    else updated_covariate_values[i] += cumulative_sum; // Add cumulative sum to existing value
	}
    }

    UNPROTECT(1);
    return updated_covariates;
}

/*
** Part 2 of the code, used for tdc
**  for each row of the baseline data (id, time1), return the row of
**  the new data (new_id, ntime2) that will provide the new data.
** Based on a last-value-carried forward rule, if the master for an id
**  had time intervals of (0,5), (5,10), (15,20) and the new data had time values
**  of  -1, 5, 11, 12, the return index would be 1, 2, and 4.  Covariates take
**  effect at the start of an interval.  Notice that obs 3 is never used.
*/

SEXP tmerge2(SEXP id2, SEXP time1x, SEXP nid2, SEXP ntime2) {
    int k; // current index
    int n_baseline, n_new;

    int *baseline_id, *new_id;
    double *baseline_time, *new_time;
    SEXP result_indices;
    int *index;

    // Get the lengths of the baseline and new data sets
    n_baseline = LENGTH(id2);  /* baseline data set */
    n_new = LENGTH(nid2); /* new data set */

    // Get pointers to the actual data in the input vectors
    baseline_id = INTEGER(id2);
    new_id = INTEGER(nid2);
    baseline_time = REAL(time1x);
    new_time = REAL(ntime2);

    // Allocate memory for the result indices (to store the index from new data for each baseline entry)
    PROTECT(result_indices = allocVector(INTSXP, n_baseline));
    index = INTEGER(result_indices);

    /*
    ** Every subject in the new data (new_id, new_time) will be found in the baseline
    ** data set (id, time1) -- if not the parent routine has already tossed
    ** them, but not every obs in id will have a representative in the new.
    ** For those we return 0, otherwise the max k such that new_id[k]== id[i]
    ** and new_time[k] <= baseline_time[i]
    **
    ** For each element in (id, time1):
    **   set index to 0
    **   walk forward in data set 2 until the newid is >= current
    **   while (id matches and newtime <= oldtime), set pointer
    **        to this row
    */
    /*
    ** Iterate over the baseline data (ids, time1), and for each entry:
    **   - Initialize the index to 0 (default when no match is found).
    **   - Traverse the new data to find the largest k such that:
    **     new_id[k] == baseline_id[i] and new_time[k] <= baseline_time[i].
    **   - Set the result index to k+1 (since R indices are 1-based).
    */
    k = 0; // Index for traversing the new data

    // Main loop over the baseline data
    for (int i = 0; i < n_baseline; i++) {
	index[i] = 0; // Default to 0, meaning no match

	// Move forward in the new data until the new ID is >= current baseline ID
	while (k < n_new && (new_id[k] < baseline_id[i]))
	{
	    k++;
	}

	// Continue traversing while new ID matches and new time is <= baseline time
        while (k < n_new && (new_id[k] == baseline_id[i]) && (new_time[k] <= baseline_time[i]))
        {
	    index[i] = k + 1; // Set index to k+1 (R's 1-based indexing)
	    k++;
        }

	// Only decrement k if it is greater than 0, to avoid going negative
	if (k > 0) k--;  /* the next obs might need the same k */
    }

    UNPROTECT(1);
    return result_indices;
}


/*
** This routine is used by surv2data to accomplish last value carried
**  forward.  It started out as a modification of tmerge2, but then
**  simplified.
** The input data has id, time, and missing yes/no.  The return value
**  is a vector of integers pointing to the replacement element, or
**  0 if there is no good replacement, e.g., the first time point for a
**  subject is missing.  
** Id is an integer.
*/
SEXP tmerge3(SEXP id2, SEXP miss2) {
    int k; // current index
    int n;
    int last_subject_id;

    int *baseline_id, *missing;
    SEXP result_indices;
    int *index;

    n = LENGTH(id2); // Number of entries in the baseline dataset
    baseline_id = INTEGER(id2);
    missing = INTEGER(miss2);  /* actually logical, but they pass as integer*/

    // Allocate memory for the result indices (to store the index from new data for each baseline entry)
    PROTECT(result_indices = allocVector(INTSXP, n));
    index = INTEGER(result_indices);

    /*
    ** The input is sorted by time within each ID.
    ** We keep track of the last non-missing index seen for each subject ID.
    ** When encountering a new ID, we reset the last valid index.
    */
    last_subject_id = -1; // Initialize to an invalid ID
    k = 0; // Last non-missing index

    // Main loop over the baseline data
    for (int i = 0; i < n; i++) {
	if (baseline_id[i] != last_subject_id) {
	    k = 0; // Reset for a new ID
	    last_subject_id = baseline_id[i];
	}
	if (missing[i] == 1) index[i] = k; // Use last valid index if missing
	else {
	    index[i] = i + 1; // Store current index (1-based)
	    k = i + 1; // Update last valid index
	}
    }

    UNPROTECT(1);
    return result_indices;
}
