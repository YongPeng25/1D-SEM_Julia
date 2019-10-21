# function myibool(ngll,nspec)
#     ibool = zeros(Int, ngll, nspec)
#     iglob = 1
#     for ispec = 1:nspec
#       for i = 1:ngll
#         if i > 1
#           iglob = iglob + 1
#         end
#         ibool[i, ispec] = iglob
#       end
#     end
#     ibool
# end



# function myibool(ngll,nspec)
    ngll = 5; nspec = 5
    ibool = zeros(Int, ngll, nspec)
    iglob = 1
    for ispec = 1:nspec
      for i = 1:ngll
        # if i > 1
        # iglob = iglob + 1
        # end
        ibool[i, ispec] = ibool[i, ispec] + 1
      end
    end
    println(ibool)
# end

# ngll = 5; nspec = 5
# ibool = myibool(ngll, nspec)
# println(ibool)
