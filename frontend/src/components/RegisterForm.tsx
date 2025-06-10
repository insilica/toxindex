import React, { useState } from "react";

const RegisterForm: React.FC = () => {
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [passwordConfirmation, setPasswordConfirmation] = useState("");
  const [error, setError] = useState("");
  const [success, setSuccess] = useState("");
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError("");
    setSuccess("");
    setLoading(true);
    try {
      const response = await fetch("/api/register", {
        method: "POST",
        headers: { "Content-Type": "application/x-www-form-urlencoded" },
        body: new URLSearchParams({
          email,
          password,
          password_confirmation: passwordConfirmation,
        }),
      });
      if (response.redirected) {
        window.location.href = response.url;
        return;
      }
      const text = await response.text();
      if (text.includes("already exists")) {
        setError("This email is already registered. Please log in.");
      } else if (text.includes("error") || text.includes("An error occurred")) {
        setError("An error occurred while creating the account. Please try again.");
      } else {
        setSuccess("Registration successful! Please check your email to verify your account.");
      }
    } catch (err) {
      setError("An error occurred. Please try again.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <form onSubmit={handleSubmit} className="bg-white p-8 rounded shadow-md w-full max-w-sm mx-auto mt-12">
      <h2 className="text-2xl font-bold mb-6 text-center">Register</h2>
      {error && <div className="text-red-600 mb-4 text-center">{error}</div>}
      {success && <div className="text-green-600 mb-4 text-center">{success}</div>}
      <input
        type="email"
        className="block w-full mb-4 p-3 border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500 text-gray-900"
        placeholder="Email"
        value={email}
        onChange={e => setEmail(e.target.value)}
        required
      />
      <input
        type="password"
        className="block w-full mb-4 p-3 border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500 text-gray-900"
        placeholder="Password"
        value={password}
        onChange={e => setPassword(e.target.value)}
        required
      />
      <input
        type="password"
        className="block w-full mb-6 p-3 border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500 text-gray-900"
        placeholder="Confirm Password"
        value={passwordConfirmation}
        onChange={e => setPasswordConfirmation(e.target.value)}
        required
      />
      <button
        type="submit"
        className="w-full bg-green-600 text-white py-3 rounded font-semibold hover:bg-green-700 transition"
        disabled={loading}
      >
        {loading ? "Registering..." : "Register"}
      </button>
      <div className="mt-4 text-center">
        <a href="/login" className="text-blue-600 hover:underline text-sm">Already have an account? Log in</a>
      </div>
    </form>
  );
};

export default RegisterForm; 