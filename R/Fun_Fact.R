#' Funfacts about birds and maybee other wildlife
#'
#' @description randomly selects a trivia from the list  of astonishing facts
#'
#' @return
#' \item{text}{a short trivia}
#'
#' @references
#'  https://www.audubon.org/news/11-fun-facts-about-owls
#'  https://www.audubon.org/news/9-fun-facts-about-turkeys
#'  https://www.audubon.org/news/9-awesome-facts-about-bird-migration
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @keywords internal
#' @export
#'
#'
#'
Fun_Fact <- function(){
    trivia <- list()
    trivia[[1]] <-
        c("Do you know this?","",
              "The probably oldest bird nest has been discovered on a cliff in Greenland.",
              "The nesting site was already used 2500 years ago and is still used by gryfalcons,",
              "the largest species of falcon in the world","Burnham et al. Ibis (2009) 151.")
    trivia[[2]] <-
        c("Do you know how well birds can see?","",
        "The Northern Hawk Owl (Surnia ulula) can detect primarily by sight a vole to eat up to a half a mile away.")
    trivia[[3]] <-
        c("Do you know in which great heights birds migrate?","",
          "Bar\u002Dheaded geese are the highest flying migratory birds,",
        "regularly reaching altitudes of up to five and a half miles above sea level",
        "while flying over the Himalayas in India.",
        "But the bird with the record for the highest altitude ever is the Ruppel\u0027s griffon vulture,",
          "which collided with a plane at 37,000 feet (that is seven miles!) in 1975",
        "and was unfortunately sucked into its jet engine." )
    trivia[[4]] <-
        c("Do you know over which distances birds migrate?","",
            "The Arctic tern has the longest migration of any bird in the world. They can fly more than 49,700 miles in a year",
          "making a round trip between their breeding grounds in the Arctic and the Antarctic, where they spend their winters.",
          " The lucky bird gets to see two summers a year! And over its lifespan of more than 30 years,",
            "the flights can add up to the equivalent of three trips to the moon and back")
    trivia[[5]] <-
        c("You know of the longest non stop flight of a bird ever recorded?","",
          "The bar\u002Dtailed godwit can fly for nearly 7,000 miles without stopping, making it the bird with the longest recorded non\u002Dstop flight.",
          "During the eight\u002Dday journey, the bird does not stop for food or rest, demonstrating jaw dropping endurance")
    x <- paste(trivia[[sample(x = 1:length(trivia),size = 1)]],collapse = " ")
    cat(stringr::str_wrap(x,width=50,exdent=2,indent=2))
}

