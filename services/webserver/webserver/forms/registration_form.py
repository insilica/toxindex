from flask import Flask
from flask_wtf import FlaskForm
from wtforms import StringField,PasswordField
from wtforms import validators
from wtforms import StringField, PasswordField, ValidationError
from wtforms.validators import DataRequired, Length, Email, EqualTo

class RegistrationForm(FlaskForm):
    email = StringField('Email', [
        DataRequired("Please enter your email address."),
        Length(min=6, max=35),
        Email("Please enter a valid email address.")
    ])
    password = PasswordField('Password', [
        DataRequired("Please enter a password.")
    ])
    password_confirmation = PasswordField('Confirm Password', [
        DataRequired("Please confirm your password."),
        EqualTo('password', message="Passwords must match.")
    ])

class UpdatePasswordForm(FlaskForm):
  password = PasswordField('password', [validators.InputRequired()])
  
class ForgotPasswordForm(FlaskForm):
  email    = StringField('my email', [validators.Length(min=6, max=35)])