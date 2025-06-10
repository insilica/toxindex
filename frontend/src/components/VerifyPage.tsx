import React from "react";

const VerifyPage: React.FC = () => {
  // Optionally, you could parse the email from the query string here
  return (
    <div className="max-w-md mx-auto mt-12 p-8 bg-white rounded shadow">
      <h2 className="text-2xl font-bold mb-4 text-center">Verify Your Email</h2>
      <p className="text-center">
        Please check your email and click the verification link to activate your account.
      </p>
    </div>
  );
};

export default VerifyPage; 