import React, { useState } from "react";

const LoginForm: React.FC = () => {
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError("");
    console.log('[LoginForm] Form submitted');
    console.log('[LoginForm] Email:', email);
    console.log('[LoginForm] Password length:', password.length);
    try {
      const response = await fetch("/api/login", {
        method: "POST",
        headers: { "Content-Type": "application/x-www-form-urlencoded" },
        body: new URLSearchParams({ email, password }),
      });
      console.log('[LoginForm] Response status:', response.status);
      if (response.redirected) {
        console.log('[LoginForm] Redirected to:', response.url);
        window.location.href = response.url;
        return;
      }
      const text = await response.text();
      console.log('[LoginForm] Response text:', text);
      if (text.includes("bad login")) {
        setError("Invalid email or password.");
      } else {
        localStorage.setItem("auth_token", "dummy_token");
        window.location.reload();
      }
    } catch (err) {
      console.error('[LoginForm] Error:', err);
      setError("An error occurred. Please try again.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <form onSubmit={handleSubmit} className="bg-white p-8 rounded shadow-md w-full max-w-sm mx-auto mt-12">
      <h2 className="text-2xl font-bold mb-6 text-center">Login</h2>
      {error && <div className="text-red-600 mb-4 text-center">{error}</div>}
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
        className="block w-full mb-6 p-3 border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500 text-gray-900"
        placeholder="Password"
        value={password}
        onChange={e => setPassword(e.target.value)}
        required
      />
      <button
        type="submit"
        className="w-full bg-green-600 text-white py-3 rounded font-semibold hover:bg-green-700 transition"
        disabled={loading}
      >
        {loading ? "Logging in..." : "Login"}
      </button>
      <div className="mt-4 text-center">
        <a href="/forgot_password" className="text-blue-600 hover:underline text-sm">Forgot password?</a>
      </div>
      <div className="mt-2 text-center">
        <a href="/register" className="text-blue-600 hover:underline text-sm">Need to create an account?</a>
      </div>
    </form>
  );
};

export default LoginForm; 