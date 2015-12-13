#' Email Notifications
#'
#' Emails a message from R when called, useful for code with long computation time
#'
#' @param subject character string for subject of email. Defaults to "Job Finished"
#' @param body optional character string for body of email. Defaults to nothing
#' @param email character string for email address to send to/from. Defaults to "jeradhoyATgmail.com"
#'
#' @return None
#'
#' @examples
#' notifyMe() 
#' notifyMe("Routing Finished")
#' notifyMe("Routing Finished", "Finished routing historical dayment data", emailaddressATexample.com)
#' notifyMe(paste("Routing Finished at", Sys.time()))
#' @name notifyMe
#' @export

notifyMe <- function(subject= "Job Finished", body=" ", email="jeradhoy@gmail.com"){
    send.mail(from = email,
	      to = email,
	      subject = paste("R -", subject),
	      body = paste(body, Sys.time()),
	      smtp = list(host.name = "aspmx.l.google.com", port = 25),
	      authenticate = FALSE,	 
	      send = TRUE)
}
