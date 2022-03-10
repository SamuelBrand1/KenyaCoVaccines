"""
        function get_rescaled_contact_matrices(filename::AbstractString)

Get three age-group contact matrices: contacts at home, work and social contacts, and, contacts
at school. The `JLD2` file at `filename` must contain the underlying contact matrices.
"""
function get_rescaled_contact_matrices(filename::AbstractString)
        M_Data = FileIO.load(filename)
        M_Kenya = M_Data["M_Kenya"]
        M_Kenya_ho = M_Data["M_Kenya_ho"]
        M_Kenya_other = M_Data["M_Kenya_other"]
        M_Kenya_work = M_Data["M_Kenya_work"]
        M_Kenya_school = M_Data["M_Kenya_school"]
        Ronescaler = 1/Real(eigvals(M_Kenya)[end])
        M_Kenya_ho .= M_Kenya_ho.*Ronescaler
        M_Kenya_other .= M_Kenya_other.*Ronescaler
        M_Kenya_work .= M_Kenya_work.*Ronescaler
        M_Kenya_school .= M_Kenya_school.*Ronescaler
        return M_Kenya_ho,M_Kenya_other,M_Kenya_school,M_Kenya_work
end

"""
        function get_population_size_matrix(filename)

Get a `NamedArray` data structure for numbers in each age group and each county of Kenya, as per 2019 census.
"""
function get_population_size_matrix(filename)
        df_pop = DataFrame(CSV.File(filename,stringtype=String))
        #df_pop = DataFrame(CSV.File(filename))
        N = NamedArray([Float64(df_pop[i, a+1][1]) for a = 1:17, i = 1:length(df_pop.county)])
        setnames!(N, vcat([string((a - 1) * 5) * "-" * string(a * 5 - 1) for a = 1:16], ["80+"]), 1)
        setnames!(N, df_pop.county, 2)
        setdimnames!(N, ["Age", "County"])
        return N
end

"""
 Function to reduce the age catorgories from 17 to 6
 This function collapses age groups into [0-19], [20-49], [50-59], [60-69], [70-79], 80+
"""

function reduce_age_categories(M_age,N_kenya,county)
    N_county = N_kenya[:,county]
    P_b =  N_county/sum(N_county)
    p_b_prime = vcat(sum(P_b[1:4]),sum(P_b[5:10]),sum(P_b[11:12]),sum(P_b[13:14]),sum(P_b[15:16]),P_b[17])
    p_vec = [[1,2,3,4],[5,6,7,8,9,10],[11,12],[13,14],[15,16],[17]]
    M_age_prime = zeros(6,6)

    for a_prime in 1:6, b_prime in 1:6
        M_age_prime[a_prime,b_prime] = sum([M_age[a,b]*P_b[b]/p_b_prime[b_prime] for a in p_vec[a_prime] for b in p_vec[b_prime]])
    end
    return  M_age_prime
end
