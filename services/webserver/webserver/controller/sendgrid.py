import sendgrid
import os
from sendgrid.helpers.mail import Mail, Email, To, Content


def send(email,subject,content):
    sg = sendgrid.SendGridAPIClient(api_key=os.getenv('SENDGRID_API_KEY'))

    from_email = Email("info@insilica.co") # TODO create toxindex email

    mail = Mail(from_email=from_email, to_emails=email, subject=subject, html_content=content)

    response = sg.send(mail)
    return response